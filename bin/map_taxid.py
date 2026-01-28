import regex as re

MAJOR_RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
mapping_ranks = {
    "acellular root": "superkingdom",
    "domain": "superkingdom"
}
forced_taxid_to_names = {
    "mycoplasma arthritidis str. 158p10p9": 243272,
    "clostridium sordelli": 1505,
    "rhodococcoides fascians": 1828
}

def get_root(taxid, desired_rank, taxdump_dict):
    """
    Walk up the taxonomy tree until desired_rank is found.

    Parameters
    ----------
    taxid : int or str
        Starting taxid
    desired_rank : str
        Rank to stop at (e.g. 'genus', 'family')
    taxdump_dict : dict
        Taxonomy dictionary keyed by taxid

    Returns
    -------
    str or None
        Taxid at the desired rank, or None if not found
    """
    taxid = str(taxid)

    while taxid in taxdump_dict:
        node = taxdump_dict[taxid]

        if node.get("rank") == desired_rank:
            return taxid

        parent = node.get("parent_taxid")

        # Stop if no parent or self-loop
        if not parent or parent == taxid:
            break

        taxid = parent

    return None



def _split_dmp_line(line):
    parts = [p.strip() for p in line.rstrip().split("|")]
    return [p for p in parts if p != ""]

def load_merged(merged_path):
    m = {}
    with open(merged_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            p = _split_dmp_line(line)
            try:
                m[int(p[0])] = int(p[1])
            except Exception:
                pass
    return m

def load_nodes(nodes_path):
    parent = {}
    rank = {}
    with open(nodes_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            p = _split_dmp_line(line)
            try:
                tid  = int(p[0])
                ptid = int(p[1])
                rnk  = p[2]
                rnk_norm = mapping_ranks.get(rnk.lower(), rnk)
                parent[tid] = ptid
                rank[tid]   = rnk_norm
            except Exception:
                pass
    return parent, rank
VALID_NAME_CLASSES = {"scientific name", "synonym", "equivalent name", "misspelling", "authority", 'basionym'}

BINOMIAL_RE = re.compile(r"^([A-Z][a-z]+)\s+([a-z][a-z0-9_-]+)")

def strip_before_paren(name: str) -> str:
    """'X (authorship...)' -> 'X' (unchanged if no '(')."""
    return re.split(r"\s*\(", name, maxsplit=1)[0].strip()

def extract_binomial(name: str) -> str | None:
    """Return 'Genus species' if name starts with a binomial, else None."""
    m = BINOMIAL_RE.match(name.strip())
    if not m:
        return None
    return f"{m.group(1)} {m.group(2)}"

def synonym_aliases(name_txt: str) -> set[str]:
    """
    Produce a set of alias strings to also store for synonyms:
      - before '('
      - binomial
    """
    aliases = set()

    # Method 1: strip parenthetical authorship
    a1 = strip_before_paren(name_txt)
    if a1:
        aliases.add(a1)

    # Method 2: extract binomial (works with or without parentheses)
    a2 = extract_binomial(name_txt)
    if a2:
        aliases.add(a2)

    # Optional: drop exact original (we already store full synonym separately)
    aliases.discard(name_txt.strip())

    return aliases


def load_taxdump(taxdump):
    taxdump_dict = {}
    with open(taxdump) as f:
        for line in f:
            # Assuming the file is tab-delimited and the columns are ordered as:
            # taxid, parent_taxid, rank, name
            try:
                parts = line.strip().split("\t")
                taxid = parts[0]
                parent_taxid = parts[2]
                rank = parts[4]
                taxdump_dict[taxid] = {
                    'parent_taxid': parent_taxid,
                    'rank': rank
                }
            except Exception as ex:
                print(f"Error parsing line: {line}, {ex}")
    f.close()
    return taxdump_dict

def load_names(names_dmp_path):
    """
    Parse names.dmp (NCBI) and return a dict mapping taxid -> scientific name.
    Only keeps the 'scientific name' class.
    """
    names_map = {}
    with open(names_dmp_path, 'r') as fh:
        for line in fh:
            # Typical format: "<taxid>\t|\t<name>\t|\t<unique name>\t|\t<name class>\t|"
            parts = [p.strip() for p in line.split('|')]
            if len(parts) < 4:
                continue
            taxid = parts[0]
            name = parts[1]
            name_class = parts[3]
            if name_class == "scientific name":
                names_map[taxid] = name
    return names_map

def load_name_records(names_path):
    """
    Return a single dict:
      records[name_lower] = {
          "taxid": int,
          "primary": bool,          # True iff scientific name
          "name_class": str,        # from names.dmp
          "name": str               # canonical/primary scientific name (filled in pass 2)
      }

    Pass 1:
      - store all names (limited to VALID_NAME_CLASSES)
      - store primary scientific name per taxid
    Pass 2:
      - for non-primary entries, set records[name_lower]["name"] to the primary scientific name for that taxid
      - if no primary found, leave "name" as the raw name (fallback)
    """
    records = {}          # name_lower -> record
    primary_by_taxid = {} # taxid -> primary scientific name (lowercased)

    with open(names_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            p = _split_dmp_line(line)
            try:
                tid = int(p[0])
                name_txt = p[1].strip()
                name_class = p[-1].strip()
            except Exception:
                continue

            if name_class not in VALID_NAME_CLASSES:
                continue

            key = name_txt.lower()
            is_primary = (name_class == "scientific name")

            # Store the full name
            if key not in records:
                records[key] = {
                    "taxid": tid,
                    "primary": is_primary,
                    "name_class": name_class,
                    "name": key
                }

            # If this is a synonym, also store the stripped base name
            if name_class in ['authority', "synonym", "basionym"]:
                for alias in synonym_aliases(name_txt):
                    alias_key = alias.lower()
                    if alias_key not in records:
                        records[alias_key] = {
                            "taxid": tid,
                            "primary": False,
                            "name_class": "synonym",
                            "name": alias_key,  # pass 2 will rewrite to primary scientific name
                        }

            # Track the primary scientific name for this taxid (first wins)
            if is_primary and tid not in primary_by_taxid:
                primary_by_taxid[tid] = key

    # Pass 2: assign canonical primary scientific name to each alt name
    for key, rec in records.items():
        tid = rec["taxid"]
        primary_name = primary_by_taxid.get(tid)

        if primary_name:
            # For primary entries, this is just itself; for alts it becomes the canonical name
            rec["name"] = primary_name
        else:
            # Fallback: no scientific name seen for this taxid in the file (rare)
            # keep whatever we already had
            rec["name"] = rec["name"]

    return records

def strip_authorship(name):
    """
    Extract base taxon name by removing authorship info.
    Example:
      'Clostridium difficile (Hall and O'Toole 1935) Prevot 1938'
      -> 'Clostridium difficile'
    """
    return re.split(r"\s*\(", name, maxsplit=1)[0].strip()

def resolve_current_taxid(tid, merged_map):
    seen = set(); cur = tid
    while cur in merged_map and cur not in seen:
        seen.add(cur); cur = merged_map[cur]
    return cur

def climb_lineage(tid, parent_map, rank_map):
    out = []; cur = tid; seen = set()
    while cur in parent_map and cur not in seen:
        seen.add(cur)
        out.append((cur, rank_map.get(cur, "no rank")))
        p = parent_map[cur]
        if p == cur:
            break
        cur = p
    return out

def extract_major_ranks(lineage_pairs):
    res = {r: "" for r in MAJOR_RANKS}
    for tid, rnk in lineage_pairs:
        r = mapping_ranks.get(rnk.lower() if isinstance(rnk, str) else rnk, rnk)
        if r in res and res[r] == "":
            res[r] = str(tid)
    return res

def map_names_for_ranks(rank_taxids, sci_name_map):
    out = {}
    for rnk, tid_str in rank_taxids.items():
        out[f"{rnk}_name"] = sci_name_map.get(int(tid_str), "") if tid_str else ""
    return out

