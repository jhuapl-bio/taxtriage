#!/usr/bin/env python3
"""Tests for bin/index_hopping.py (run: python3 bin/test_index_hopping.py)."""
import copy
import json
import os
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import index_hopping as ih  # noqa: E402

# The layout from issue #344 (tab-separated, no hopped rows).
USER_TSV = (
    "Lane\tSampleID\tSample_Project\tindex\tindex2\t# Reads\t% of Hopped Reads\t% of All Reads\n"
    "1\tIDR2600053924-01-01\tMTB-PHOENIX-11955\tTGCCGGTCAG\tTTGTATCAGG\t1861265\t\t0.003531\n"
    "1\tIDR2600055130-01-01\tMTB-PHOENIX-11955\tCACTCAATTC\tCACAACTTAA\t1287738\t\t0.002443\n"
    "1\tIDR2600055131-01-01\tMTB-PHOENIX-11955\tTCTCACACGC\tTCACTGTCCG\t2116690\t\t0.004016\n"
    "1\tIDR2600055132-01-01\tMTB-PHOENIX-11955\tTCAATGGAGA\tGTTGTAGTGC\t3024655\t\t0.005739\n"
)

# bcl-convert style CSV with hopped combinations (blank SampleID).
BCL_CSV = (
    "Lane,SampleID,index,index2,# Reads,% of Hopped Reads,% of All Reads\n"
    "1,A,AAAAAAAA,CCCCCCCC,1000000,,0.5\n"
    "1,B,GGGGGGGG,TTTTTTTT,1000000,,0.5\n"
    "1,,AAAAAAAA,TTTTTTTT,1500,0.75,0.00075\n"
    "1,,GGGGGGGG,CCCCCCCC,500,0.25,0.00025\n"
    "2,A,AAAAAAAA,CCCCCCCC,10,,0.5\n"
    "2,C,ACACACAC,GTGTGTGT,10,,0.5\n"
    "2,,AAAAAAAA,GTGTGTGT,1,1,0.05\n"
)


def write_xlsx(path, rows):
    """Minimal XLSX with inline strings (no openpyxl needed)."""
    def col(i):
        s = ""
        i += 1
        while i:
            i, r = divmod(i - 1, 26)
            s = chr(65 + r) + s
        return s
    xml_rows = []
    for ri, row in enumerate(rows, 1):
        cells = []
        for ci, v in enumerate(row):
            ref = f"{col(ci)}{ri}"
            if isinstance(v, (int, float)):
                cells.append(f'<c r="{ref}"><v>{v}</v></c>')
            else:
                cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{v}</t></is></c>')
        xml_rows.append(f'<row r="{ri}">{"".join(cells)}</row>')
    sheet = ('<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/'
             'spreadsheetml/2006/main"><sheetData>' + "".join(xml_rows) + "</sheetData></worksheet>")
    wb = ('<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="http://schemas.openxmlformats.org/'
          'spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/'
          'relationships"><sheets><sheet name="Index Hopping" sheetId="1" r:id="rId1"/></sheets></workbook>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/'
            'package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("[Content_Types].xml", "<Types/>")
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet)


def node(key, name, reads, members=None, tass=50.0):
    n = {"key": key, "name": name, "numreads": float(reads), "tass_score": tass,
         "passes_threshold": tass >= 50}
    if members is not None:
        n["members"] = members
    return n


def sample(name, org_reads):
    """org_reads: {strain_key: reads}; all under one species/genus."""
    strains = [node(k, f"Org {k}", r) for k, r in org_reads.items()]
    sp = node("sp1", "Species 1", sum(org_reads.values()), strains)
    gn = node("g1", "Genus 1", sum(org_reads.values()), [sp])
    return {"metadata": {"sample_name": name}, "organisms": [gn]}


class TestParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _w(self, name, text):
        p = os.path.join(self.tmp, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def test_user_tsv(self):
        hf = ih.parse_hopping_file(self._w("hop.tsv", USER_TSV))
        self.assertEqual(list(hf["lanes"]), ["1"])
        L = hf["lanes"]["1"]
        self.assertEqual(len(L["samples"]), 4)
        self.assertEqual(L["hopped"], [])
        self.assertEqual(L["samples"]["IDR2600055132-01-01"]["reads"], 3024655)
        self.assertTrue(hf["has_index2"])

    def test_same_table_as_txt_and_csv(self):
        a = ih.parse_hopping_file(self._w("hop.txt", USER_TSV))
        b = ih.parse_hopping_file(self._w("hop.csv", USER_TSV.replace("\t", ",")))
        self.assertEqual(a["lanes"], b["lanes"])

    def test_bcl_csv_hops(self):
        hf = ih.parse_hopping_file(self._w("Index_Hopping_Counts.csv", BCL_CSV))
        self.assertEqual(sorted(hf["lanes"]), ["1", "2"])
        self.assertEqual(len(hf["lanes"]["1"]["hopped"]), 2)
        lr = ih.LaneRates("1", hf["lanes"]["1"])
        self.assertAlmostEqual(lr.lane_rate, 2000 / 2002000)
        rate, src = lr.rate("A", "B", 0.001)
        self.assertEqual(src, "pair")
        self.assertAlmostEqual(rate, 2001 / 2000000)

    def test_xlsx(self):
        rows = [r.split("\t") for r in USER_TSV.strip().split("\n")]
        rows = [rows[0]] + [[int(c) if c.isdigit() else c for c in r] for r in rows[1:]]
        p = os.path.join(self.tmp, "hop.xlsx")
        write_xlsx(p, [["Index Hopping Report"], []] + rows)  # title rows above the header
        hf = ih.parse_hopping_file(p)
        tsv = ih.parse_hopping_file(self._w("hop.tsv", USER_TSV))
        self.assertEqual(hf["lanes"], tsv["lanes"])

    def test_thousands_separators(self):
        txt = USER_TSV.replace("1861265", "1,861,265")
        hf = ih.parse_hopping_file(self._w("hop.tsv", txt))
        self.assertEqual(hf["lanes"]["1"]["samples"]["IDR2600053924-01-01"]["reads"], 1861265)

    def test_bad_header(self):
        with self.assertRaises(ValueError):
            ih.parse_hopping_file(self._w("x.csv", "a,b,c\n1,2,3\n"))

    def test_poisson(self):
        self.assertAlmostEqual(ih.poisson_sf(1, 0.4), 1 - 2.718281828 ** -0.4, places=6)
        self.assertEqual(ih.poisson_sf(0, 1), 1.0)
        self.assertEqual(ih.poisson_sf(1000, 1), 0.0)


class TestAnalysis(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        p = os.path.join(self.tmp, "hop.tsv")
        with open(p, "w") as fh:
            fh.write(USER_TSV)
        self.hf = [ih.parse_hopping_file(p)]

    def run_case(self, samples, **kw):
        args = dict(default_rate=0.001, factor=2.0, alpha=0.01, log=lambda *_: None)
        args.update(kw)
        return ih.analyse(samples, self.hf, None, {}, **args)

    def test_flags_low_level_hop(self):
        s = {
            "IDR2600053924-01-01": sample("d", {"mtb": 50000}),   # donor
            "IDR2600055130-01-01": sample("r", {"mtb": 40}),      # 40 <= Pois(2*50)
            "IDR2600055131-01-01": sample("x", {"mtb": 5000}),    # real
        }
        rows = self.run_case(s)
        rec = s["IDR2600055130-01-01"]["organisms"][0]["members"][0]["members"][0]["index_hop"]
        self.assertTrue(rec["suspect"])
        self.assertEqual(rec["top_donor"], "IDR2600053924-01-01")
        self.assertEqual(rec["rate_source"], "default")
        self.assertAlmostEqual(rec["expected_reads"], 50000 * 0.001 + 5000 * 0.001)
        real = s["IDR2600055131-01-01"]["organisms"][0]["members"][0]["members"][0]["index_hop"]
        self.assertFalse(real["suspect"])
        donor = s["IDR2600053924-01-01"]["organisms"][0]["members"][0]["members"][0]["index_hop"]
        self.assertFalse(donor["suspect"])
        self.assertEqual(donor["donors"], [])
        # flagged at every level (strain, species, genus) of the recipient
        self.assertEqual(sorted(r["level"] for r in rows if r["sample"] == "IDR2600055130-01-01"),
                         ["genus", "species", "strain"])
        meta = s["IDR2600055130-01-01"]["metadata"]["index_hopping"]
        self.assertTrue(meta["applied"])
        self.assertEqual(meta["n_suspect"]["strain"], 1)

    def test_does_not_change_scores(self):
        s = {"IDR2600053924-01-01": sample("d", {"mtb": 50000}),
             "IDR2600055130-01-01": sample("r", {"mtb": 40})}
        before = copy.deepcopy(s)
        self.run_case(s)

        def strip(o):
            if isinstance(o, dict):
                return {k: strip(v) for k, v in o.items() if k not in ("index_hop", "index_hopping", "index_hop_source")}
            if isinstance(o, list):
                return [strip(v) for v in o]
            return o
        self.assertEqual(strip(s), strip(before))

    def test_source_side_note(self):
        s = {
            "IDR2600053924-01-01": sample("d", {"mtb": 50000}),   # donor: 50 expected -> noted
            "IDR2600055130-01-01": sample("r", {"mtb": 40}),
            "IDR2600055131-01-01": sample("x", {"mtb": 300}),     # 0.3 expected -> below MIN_SOURCE_READS
        }
        self.run_case(s)
        d = s["IDR2600053924-01-01"]["organisms"][0]["members"][0]["members"][0]
        self.assertEqual(d["index_hop_source"]["recipients"][0]["sample"], "IDR2600055130-01-01")
        self.assertAlmostEqual(d["index_hop_source"]["expected_out"], 50.0)
        self.assertNotIn("suspect", d["index_hop_source"])
        x = s["IDR2600055131-01-01"]["organisms"][0]["members"][0]["members"][0]
        self.assertNotIn("index_hop_source", x)
        self.assertEqual(s["IDR2600053924-01-01"]["metadata"]["index_hopping"]["n_source"]["strain"], 1)

    def test_unmatched_sample(self):
        s = {"IDR2600053924-01-01": sample("d", {"mtb": 50000}),
             "Other_Run_Sample": sample("o", {"mtb": 3})}
        self.run_case(s)
        self.assertFalse(s["Other_Run_Sample"]["metadata"]["index_hopping"]["applied"])
        self.assertNotIn("index_hop", s["Other_Run_Sample"]["organisms"][0])

    def test_bcl_suffix_and_case(self):
        s = {"idr2600053924_01_01_S1": sample("d", {"mtb": 50000}),
             "IDR2600055130-01-01": sample("r", {"mtb": 10})}
        self.run_case(s)
        self.assertTrue(s["idr2600053924_01_01_S1"]["metadata"]["index_hopping"]["applied"])
        self.assertTrue(s["IDR2600055130-01-01"]["organisms"][0]["members"][0]["members"][0]
                        ["index_hop"]["suspect"])

    def test_override_sample_id(self):
        s = {"donor": sample("d", {"mtb": 50000}), "recip": sample("r", {"mtb": 10})}
        ov = {"donor": "IDR2600053924-01-01", "recip": "IDR2600055130-01-01"}
        ih.analyse(s, self.hf, None, ov, 0.001, 2.0, 0.01, log=lambda *_: None)
        self.assertTrue(s["recip"]["organisms"][0]["members"][0]["members"][0]["index_hop"]["suspect"])

    def test_separate_files_do_not_mix(self):
        p2 = os.path.join(self.tmp, "hop2.tsv")
        with open(p2, "w") as fh:
            fh.write(USER_TSV.replace("IDR2600055130-01-01", "OTHER-01"))
        hfs = self.hf + [ih.parse_hopping_file(p2)]
        s = {"IDR2600053924-01-01": sample("d", {"mtb": 50000}), "OTHER-01": sample("r", {"mtb": 10})}
        ih.analyse(s, hfs, {"IDR2600053924-01-01": [0], "OTHER-01": [1]}, {}, 0.001, 2.0, 0.01,
                   log=lambda *_: None)
        rec = s["OTHER-01"]["organisms"][0]["members"][0]["members"][0]["index_hop"]
        self.assertFalse(rec["suspect"])
        self.assertEqual(rec["donors"], [])

    def test_pair_rate_and_lanes(self):
        p = os.path.join(self.tmp, "bcl.csv")
        with open(p, "w") as fh:
            fh.write(BCL_CSV)
        hf = [ih.parse_hopping_file(p)]
        s = {"A": sample("A", {"o": 100000}), "B": sample("B", {"o": 150}), "C": sample("C", {"o": 150})}
        ih.analyse(s, hf, None, {}, 0.001, 2.0, 0.01, log=lambda *_: None)
        b = s["B"]["organisms"][0]["members"][0]["members"][0]["index_hop"]
        c = s["C"]["organisms"][0]["members"][0]["members"][0]["index_hop"]
        self.assertEqual(b["rate_source"], "pair")         # lane 1: A/B hopped rows
        self.assertTrue(b["suspect"])                       # 100000 * ~1e-3 = ~100 expected
        self.assertEqual(c["rate_source"], "pair")         # lane 2: A/C pair row (1 read)
        self.assertAlmostEqual(c["top_donor_rate"], 2 / 20)  # (1 + 1) / (10 + 10)

    def test_cli_end_to_end(self):
        d = os.path.join(self.tmp, "in")
        os.makedirs(d)
        paths = []
        for sid, reads in (("IDR2600053924-01-01", 50000), ("IDR2600055130-01-01", 20)):
            p = os.path.join(d, f"{sid}.paths.json")
            with open(p, "w") as fh:
                json.dump(sample(sid, {"mtb": reads}), fh)
            paths.append(p)
        smap = os.path.join(self.tmp, "map.json")
        with open(smap, "w") as fh:
            json.dump({"samples": {"IDR2600053924-01-01": {"files": [0]},
                                   "IDR2600055130-01-01": {"files": [0]}}}, fh)
        out = os.path.join(self.tmp, "out")
        tsv = os.path.join(self.tmp, "s.tsv")
        rep = os.path.join(self.tmp, "r.json")
        ih.main(["-i", *paths, "--hopping", self.hf[0]["path"], "--sample-map", smap,
                 "-o", out, "--summary", tsv, "--report", rep])
        self.assertEqual(sorted(os.listdir(out)), sorted(os.path.basename(p) for p in paths))
        with open(tsv) as fh:
            lines = fh.read().strip().split("\n")
        self.assertEqual(len(lines), 4)  # header + strain/species/genus
        with open(rep) as fh:
            self.assertIn("IDR2600055130-01-01", json.load(fh)["samples"])


if __name__ == "__main__":
    unittest.main(verbosity=1)
