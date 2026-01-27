import re

_GT_RE = re.compile(r"^(.*)_\d+_\d+$")

def ground_truth_from_read_id(read_id: str) -> str:
    """
    Extract ground-truth reference accession from read_id like:
      NZ_CP076401.1_27_5  -> NZ_CP076401.1
    Falls back to read_id if pattern doesn't match.
    """
    m = _GT_RE.match(read_id)
    return m.group(1) if m else read_id

def normalize_ref_id(ref: str) -> str:
    return (ref or "").strip()

def strip_version(ref: str) -> str:
    # NZ_CP076401.1 -> NZ_CP076401
    r = normalize_ref_id(ref)
    return r.split(".", 1)[0] if "." in r else r
import pysam
from collections import defaultdict
def build_ref_metadata(reference_hits: dict):
    meta = {}
    for ref, v in reference_hits.items():
        org = v.get("name") or v.get("assemblyname") or ref
        tax = str(v.get("taxid")) if v.get("taxid") is not None else None

        r0 = normalize_ref_id(ref)
        meta[r0] = {"taxid": tax, "organism": org}

        r1 = strip_version(r0)
        if r1 and r1 not in meta:
            meta[r1] = {"taxid": tax, "organism": org}

    return meta
def collect_read_alignments(bam_path: str, alignments_to_remove=None):
    """
    Returns:
      read_to_refs: dict(read_id -> set(refs aligned))
      read_to_seq:  dict(read_id -> seq)
      all_reads:    set(read_id) includes unmapped too
    alignments_to_remove format:
      alignments_to_remove[read_id] = set(refs_to_remove)   (or something equivalent)
    """
    read_to_refs = defaultdict(set)
    read_to_seq  = {}
    all_reads    = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # until_eof=True ensures we see unmapped reads too
        for read in bam.fetch(until_eof=True):
            rid = read.query_name
            all_reads.add(rid)

            if rid not in read_to_seq and read.query_sequence:
                read_to_seq[rid] = read.query_sequence

            if read.is_unmapped:
                continue

            ref = read.reference_name

            # Apply logical filtering if we don't have an actual filtered BAM
            if alignments_to_remove:
                # support either:
                #   alignments_to_remove[rid] is a set of refs
                # or alignments_to_remove[rid] contains ref keys etc.
                rm = alignments_to_remove.get(rid) or alignments_to_remove.get(ref)  # be conservative
                if isinstance(rm, set) and ref in rm:
                    continue
                if isinstance(alignments_to_remove.get(rid, None), (set, dict)) and ref in alignments_to_remove.get(rid, set()):
                    continue

            read_to_refs[rid].add(ref)

    return read_to_refs, read_to_seq, all_reads
import pandas as pd
def build_confusion_dataframe(read_to_refs, all_reads, ref_meta=None, prefix="NEW_"):
    read_to_truth = {rid: ground_truth_from_read_id(rid) for rid in all_reads}

    refs = set(read_to_truth.values())
    for rid, rset in read_to_refs.items():
        refs.update(rset)

    total = len(all_reads)
    rows = []

    for ref in sorted(refs):
        tp = fp = fn = 0

        for rid in all_reads:
            truth = read_to_truth[rid]
            aligned = ref in read_to_refs.get(rid, set())

            if truth == ref and aligned:
                tp += 1
            elif truth != ref and aligned:
                fp += 1
            elif truth == ref and not aligned:
                fn += 1

        tn = total - tp - fp - fn
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall    = tp / (tp + fn) if (tp + fn) else 0.0
        f1        = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

        md = (ref_meta or {}).get(ref, {})
        rows.append({
            "Reference": ref,
            "Reference_TaxID": md.get("taxid"),
            "Reference_Organism": md.get("organism"),
            f"{prefix}TP": tp, f"{prefix}FP": fp, f"{prefix}FN": fn, f"{prefix}TN": tn,
            f"{prefix}Precision": precision,
            f"{prefix}Recall": recall,
            f"{prefix}F1": f1,
            f"{prefix}TruthSupport": sum(1 for rid in all_reads if read_to_truth[rid] == ref),
            f"{prefix}PredSupport": sum(1 for rid in all_reads if ref in read_to_refs.get(rid, set())),
        })

    return pd.DataFrame(rows)

def join_meta(refs, ref_meta):
    refs_sorted = sorted(refs)
    taxids = []
    orgs = []
    for r in refs_sorted:
        t, o = meta_for_ref(r, ref_meta)
        taxids.append("" if t is None else str(t))
        orgs.append("" if o is None else str(o))
    return ",".join(taxids), ",".join(orgs)

def join_taxids(refs, ref_meta=None):
    ref_meta = ref_meta or {}
    out = []
    for r in sorted(refs):
        md = ref_meta.get(r) or ref_meta.get(r.split(".", 1)[0]) or {}
        t = md.get("taxid")
        if t not in [None, "", "None"]:
            out.append(str(t))
    return ",".join(out)

def join_orgs(refs, ref_meta=None):
    ref_meta = ref_meta or {}
    out = []
    for r in sorted(refs):
        md = ref_meta.get(r) or ref_meta.get(r.split(".", 1)[0]) or {}
        o = md.get("organism")
        if o not in [None, "", "None"]:
            out.append(str(o))
    return ",".join(out)

def meta_for_ref(ref, ref_meta=None):
    ref_meta = ref_meta or {}
    if not ref:
        return None, None
    md = ref_meta.get(ref) or ref_meta.get(ref.split(".", 1)[0]) or {}
    return md.get("taxid"), md.get("organism")

def build_false_positive_reads_df(read_to_refs_original, read_to_refs_filtered, read_to_seq, all_reads, ref_meta=None):
    rows = []
    for rid in sorted(all_reads):
        truth = ground_truth_from_read_id(rid)
        truth_taxid, truth_org = meta_for_ref(truth, ref_meta)

        refs_f = read_to_refs_filtered.get(rid, set())
        refs_o = read_to_refs_original.get(rid, set())

        # "still FP" if any remaining aligned ref != truth
        non_truth_refs = sorted([r for r in refs_f if r != truth])
        if not non_truth_refs:
            continue

        rows.append({
            "read_id": rid,

            "truth_ref": truth,
            "truth_taxid": truth_taxid,
            "truth_organism": truth_org,

            # OLD (original) — keep refs + taxids + orgs
            "aligned_refs_original": ",".join(sorted(refs_o)) if refs_o else "",
            "aligned_taxids_original": join_taxids(refs_o, ref_meta) if refs_o else "",
            "aligned_orgs_original": join_orgs(refs_o, ref_meta) if refs_o else "",
            "n_refs_original": len(refs_o),

            # NEW (filtered) — keep refs + taxids (no orgs)
            "aligned_refs_filtered": ",".join(sorted(refs_f)) if refs_f else "",
            "aligned_taxids_filtered": join_taxids(refs_f, ref_meta) if refs_f else "",
            "n_refs_filtered": len(refs_f),

            # Remaining FP refs after filtering (+ optional taxids)
            "remaining_fp_refs": ",".join(non_truth_refs),
            "remaining_fp_taxids": join_taxids(non_truth_refs, ref_meta) if non_truth_refs else "",
            "n_remaining_fp_refs": len(non_truth_refs),

            "sequence": read_to_seq.get(rid, ""),
        })

    return pd.DataFrame(rows)


def write_confusion_xlsx(
    out_xlsx_path: str,
    confusion_df: pd.DataFrame,
    fp_reads_df: pd.DataFrame
):
    with pd.ExcelWriter(out_xlsx_path, engine="openpyxl") as writer:
        confusion_df.to_excel(writer, sheet_name="confusion_matrix", index=False)
        fp_reads_df.to_excel(writer, sheet_name="remaining_false_positives", index=False)
