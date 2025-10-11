"""
A module providing utilities for genomic data transformation, enrichment analysis,
and shared enrichment scoring.

This module includes functionality for mapping protein identifiers to gene
identifiers, performing statistical over-representation analyses, term
database loaders, and calculating Shared Enrichment Scores (SES).
"""
from typing import Dict, List, Set, Tuple, Optional, Union
from pathlib import Path
import pandas as pd
import numpy as np
from math import comb
from typing import Sequence
import sqlite3
import requests
from io import StringIO
import argparse
import Utils
# Default SQLite path (change if your DB lives elsewhere)
db_path = Path("Data/DB.db")

def _normalize_ensembl_id(x: str) -> str:
    """Strip version suffix (e.g., ENSP00000354587.2 -> ENSP00000354587)."""
    if not isinstance(x, str) or not x:
        return ""
    x = x.strip()
    # keep only the part before the first dot
    core = x.split(".")[0]
    return core


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (name,))
    return cur.fetchone() is not None

###########################################
# 1) ID mapping: ENSP -> ENSG
###########################################
def ensp_list_to_ensg(ensps: Sequence[str],
                      db_path_local: str | Path = None,
                      distinct: bool = True,
                      with_meta: bool = False) -> List[str] | pd.DataFrame:
    """
    Batch mapping in the reverse direction: map a list of ENSP IDs to ENSG IDs.
    - Normalizes IDs (drops version suffix)
    - Chunks queries to respect SQLite variable limits

    Parameters
    ----------
    ensps : sequence of str
        Input protein IDs (ENSP...).
    db_path_local : str | Path | None
        SQLite DB path. Defaults to module-level db_path.
    distinct : bool
        Drop duplicate ENSG per protein when True.
    with_meta : bool
        If True, return a DataFrame with columns [ensg, ensp, hgnc_symbol, uniprot].
        If False, return a flat list of ENSG strings.
    """
    if db_path_local is None:
        db_path_local = db_path
    if not ensps:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    # Normalize and unique the inputs
    ensp_norm = []
    seen = set()
    for p in ensps:
        n = _normalize_ensembl_id(p)
        if n and n not in seen:
            seen.add(n)
            ensp_norm.append(n)
    if not ensp_norm:
        return [] if not with_meta else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])

    # Chunk to avoid SQLite parameter limits (default ~999)
    chunk_size = 800
    frames = []
    conn = sqlite3.connect(db_path_local)
    try:
        # If mapping table/view is missing, return ENSP passthrough (no crash)
        if not _table_exists(conn, "Ensembl_xref"):
            if with_meta:
                return pd.DataFrame({
                    "ensg": [],
                    "ensp": ensp_norm,
                    "hgnc_symbol": [],
                    "uniprot": []
                })
            return ensp_norm

        for i in range(0, len(ensp_norm), chunk_size):
            chunk = ensp_norm[i:i+chunk_size]
            placeholders = ",".join(["?"] * len(chunk))
            sql = f'''
                SELECT ensg, ensp, hgnc_symbol, uniprot
                FROM Ensembl_xref
                WHERE ensp IN ({placeholders}) AND ensg IS NOT NULL
            '''
            frames.append(pd.read_sql_query(sql, conn, params=chunk))
    finally:
        conn.close()

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["ensg","ensp","hgnc_symbol","uniprot"])
    if df.empty:
        return [] if not with_meta else df
    if distinct:
        df = df.drop_duplicates(subset=["ensg","ensp"]).reset_index(drop=True)
    if with_meta:
        return df
    return df["ensg"].dropna().astype(str).unique().tolist()


###########################################
# 2) Over-representation analysis (ORA)
###########################################
# NOTE: Legacy offline ORA utilities (GMT-based). Kept for reproducibility, not used in the default main flow.

def bh_fdr(pvals: List[float]) -> List[float]:
    """Benjamini–Hochberg FDR correction on a 1-D list of p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n+1)
    q = p * n / ranks
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty_like(q_sorted)
    out[order] = q_sorted
    return np.clip(out, 0, 1).tolist()

def hypergeom_right_tail(k: int, M: int, n: int, N: int) -> float:
    """
    Right-tail hypergeometric p-value: P[X >= k] for X~Hypergeom(N, M, n).
    N: population size, M: successes in population, n: draws, k: observed overlap.
    """
    if k > n or k > M:
        return 0.0
    num = 0
    for i in range(k, min(M, n) + 1):
        num += comb(M, i) * comb(N - M, n - i)
    den = comb(N, n)
    return float(num) / float(den) if den > 0 else 1.0

def enrich_one(genes_ensg: Set[str],
               term2genes: Dict[str, Set[str]],
               universe_size: Optional[int] = None) -> pd.DataFrame:
    """
    Compute ORA table for one ENSG set against term2genes (terms expressed in ENSG).
    Returns DataFrame with columns: ['term','k','M','n','N','p','q'] sorted by (q, p, -k)
    """
    # background as union of all term genes and the query, unless universe_size is provided
    bg = set().union(*term2genes.values()) | set(genes_ensg)
    N = len(bg) if universe_size is None else int(universe_size)
    n = len(genes_ensg)
    rows = []
    for term, tgenes in term2genes.items():
        M = len(tgenes)
        k = len(genes_ensg & tgenes)
        p = hypergeom_right_tail(k, M, n, N)
        rows.append({'term': term, 'k': k, 'M': M, 'n': n, 'N': N, 'p': p})
    qvals = bh_fdr([r['p'] for r in rows])
    for r, q in zip(rows, qvals):
        r['q'] = q
    df = pd.DataFrame(rows).sort_values(['q','p','k'], ascending=[True, True, False])
    return df

def top_terms(enrich_df: pd.DataFrame, top_n: int = 20, q_cut: float = 0.05) -> List[str]:
    """Pick top-N term IDs by q-value (≤ q_cut)."""
    sub = enrich_df[enrich_df['q'] <= q_cut].copy()
    return sub['term'].head(top_n).tolist()

###########################################
# 2b) Term database loaders (GMT, TSV)
###########################################

def load_gmt(gmt_path: str | Path,
             *,
             id_kind: str = "ENSG",
             drop_empty: bool = True) -> Dict[str, Set[str]]:
    """
    Read a GMT file (tab-separated): term \t desc \t gene1 \t gene2 ...
    Returns {term_id: {gene_ids...}}. Gene IDs are normalized and filtered to `id_kind` prefix when possible.
    """
    term2genes: Dict[str, Set[str]] = {}
    with open(gmt_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = set(_normalize_ensembl_id(g) for g in parts[2:])
            if id_kind:
                genes = {g for g in genes if g.startswith(id_kind)}
            if genes or not drop_empty:
                term2genes[term] = genes
    return term2genes


def load_term2genes(source: Union[str, Path, Dict[str, Set[str]]]) -> Dict[str, Set[str]]:
    """
    Convenience: if `source` is a path-like -> read GMT; if already a dict, pass-through.
    """
    if isinstance(source, (str, Path)):
        return load_gmt(source)
    return dict(source)

###########################################
# 3) Shared Enrichment Score (SES)
###########################################

def ses_from_ranked_terms(ranked_a: List[str], ranked_b: List[str], N: int = None) -> float:
    """
    SES = average cumulative Jaccard across k=1..N prefixes of two ranked lists.
    """
    if N is None:
        N = min(len(ranked_a), len(ranked_b))
    N = max(1, N)
    seen_a, seen_b = set(), set()
    vals = []
    for i in range(N):
        seen_a.add(ranked_a[i])
        seen_b.add(ranked_b[i])
        inter = len(seen_a & seen_b)
        union = len(seen_a | seen_b)
        vals.append(inter/union if union>0 else 0.0)
    return float(np.mean(vals))

###########################################
# 4) High-level pipeline
###########################################

def ses_from_ensp_lists(
    herb_ensp: List[str],
    disease_ensp: List[str],
    term2genes: Dict[str, Set[str]],
    top_n: int = 20,
    q_cut: float = 0.05,
    universe_size: int = None,
    db_path_local: Optional[str] = None,
    distinct: bool = True
) -> dict:
    """
    Full pipeline (DB-backed):
      ENSP lists -> ENSG (via SQLite mapping) -> ORA -> top terms -> SES
    Returns dict with SES and intermediate tables.
    """
    herb_ensg = set(ensp_list_to_ensg(herb_ensp, db_path_local=db_path_local, distinct=distinct))
    disease_ensg = set(ensp_list_to_ensg(disease_ensp, db_path_local=db_path_local, distinct=distinct))

    enrich_herb = enrich_one(herb_ensg, term2genes, universe_size=universe_size)
    enrich_dis  = enrich_one(disease_ensg, term2genes, universe_size=universe_size)

    ranked_herb = top_terms(enrich_herb, top_n=top_n, q_cut=q_cut)
    ranked_dis  = top_terms(enrich_dis,  top_n=top_n, q_cut=q_cut)

    N = min(len(ranked_herb), len(ranked_dis))
    ses = ses_from_ranked_terms(ranked_herb, ranked_dis, N=N) if N>0 else 0.0

    return {
        'SES': ses,
        'ranked_terms_herb': ranked_herb,
        'ranked_terms_disease': ranked_dis,
        'enrich_herb': enrich_herb,
        'enrich_disease': enrich_dis
    }

###########################################
# 5) Cytoscape-like entry point (ENSP in)
###########################################

def cytoscape_enrich_from_ensps(
    *,
    query_ensps_a: List[str],
    query_ensps_b: Optional[List[str]] = None,
    term_source: Union[str, Path, Dict[str, Set[str]]],
    db_path_local: Optional[str | Path] = None,
    top_n: int = 20,
    q_cut: float = 0.05,
    universe_size: Optional[int] = None,
    distinct: bool = True
) -> Dict[str, Union[float, List[str], pd.DataFrame]]:
    """
    Minimal, Cytoscape-like workflow.
    - Input: ENSP codes (like STRING/Cytoscape).
    - Step1: map ENSP -> ENSG via SQLite `Ensembl_xref(ensg, ensp, hgnc_symbol, uniprot)`.
    - Step2: ORA vs provided term database (GMT or dict).
    - Step3: If `query_ensps_b` given, compute Shared Enrichment Score (SES) between A and B.

    Returns dict with keys: 'enrich_a', ['enrich_b'], ['ranked_terms_a','ranked_terms_b'], and optional 'SES'.
    """
    term2genes = load_term2genes(term_source)

    a_ensg = set(ensp_list_to_ensg(query_ensps_a, db_path_local=db_path_local or db_path, distinct=distinct))
    enrich_a = enrich_one(a_ensg, term2genes, universe_size=universe_size)
    ranked_a = top_terms(enrich_a, top_n=top_n, q_cut=q_cut)

    out: Dict[str, Union[float, List[str], pd.DataFrame]] = {
        'enrich_a': enrich_a,
        'ranked_terms_a': ranked_a,
    }

    if query_ensps_b is not None:
        b_ensg = set(ensp_list_to_ensg(query_ensps_b, db_path_local=db_path_local or db_path, distinct=distinct))
        enrich_b = enrich_one(b_ensg, term2genes, universe_size=universe_size)
        ranked_b = top_terms(enrich_b, top_n=top_n, q_cut=q_cut)
        N = min(len(ranked_a), len(ranked_b))
        ses = ses_from_ranked_terms(ranked_a, ranked_b, N=N) if N > 0 else 0.0
        out.update({
            'enrich_b': enrich_b,
            'ranked_terms_b': ranked_b,
            'SES': ses
        })

    return out


###########################################
# 6) STRING API enrichment (ENSP in) + SES
###########################################


STRING_API_BASE = "https://string-db.org/api"

# Common functional categories across STRING enrichment
FUNCTIONAL_CATEGORIES: Set[str] = {"Process", "Function", "Component", "KEGG", "Reactome", "WikiPathways"}


def _string_join_identifiers(ids: Sequence[str]) -> str:
    """STRING API prefers %0d (CR) as the join delimiter for identifiers."""
    return "%0d".join(ids)


def string_get_string_ids(identifiers: Sequence[str], species: int = 9606, *, caller: str = "local_script", session: Optional[requests.Session] = None) -> Dict[str, str]:
    """
    Map input IDs (e.g., ENSP...) to STRING IDs (e.g., 9606.ENSP00000...).
    Returns {input_id: stringId}. If an ID fails to map, it is omitted.
    """
    if not identifiers:
        return {}
    url = f"{STRING_API_BASE}/json/get_string_ids"
    params = {
        "identifiers": _string_join_identifiers(list(identifiers)),
        "species": species,
        "caller_identity": caller,
    }
    own_session = False
    s = session
    if s is None:
        s = requests.Session()
        own_session = True
    try:
        r = s.post(url, data=params, timeout=30)
    finally:
        if own_session:
            s.close()
    r.raise_for_status()
    out: Dict[str, str] = {}
    for rec in r.json():
        # For each input, STRING may return multiple matches; pick highest score (first)
        inp = rec.get("queryItem") or rec.get("input_string", "")
        sid = rec.get("stringId")
        if inp and sid and inp not in out:
            out[inp] = sid
    return out


def string_enrichment_from_ids(identifiers: Sequence[str], species: int = 9606, *, ensure_string_ids: bool = True, caller: str = "local_script", output_format: str = "tsv", session: Optional[requests.Session] = None) -> pd.DataFrame:
    """
    Run STRING functional enrichment for a list of protein IDs.
    - If ensure_string_ids=True, first map to STRING IDs via get_string_ids.
    - Returns a DataFrame of enrichment results (GO/KEGG/Reactome/domains...).
    """
    ids = list(identifiers)
    if ensure_string_ids:
        mapping = string_get_string_ids(ids, species=species, caller=caller, session=session)
        ids = list(mapping.values())
        if not ids:
            return pd.DataFrame()

    method = "enrichment"
    url = f"{STRING_API_BASE}/{output_format}/{method}"
    params = {
        "identifiers": _string_join_identifiers(ids),
        "species": species,
        "caller_identity": caller,
    }
    own_session = False
    s = session
    if s is None:
        s = requests.Session()
        own_session = True
    try:
        r = s.post(url, data=params, timeout=60)
    finally:
        if own_session:
            s.close()
    r.raise_for_status()

    if output_format == "json":
        return pd.json_normalize(r.json())
    # default: tsv
    return pd.read_csv(StringIO(r.text), sep="\t")


def _string_term_key(df_row: pd.Series) -> str:
    """Create a stable term key to avoid collisions across categories."""
    # Columns present in STRING enrichment TSV: 'category', 'term', 'term_id', ...
    cat = str(df_row.get("category", "")).strip()
    tid = str(df_row.get("term_id", "")).strip()
    if tid:
        return f"{cat}:{tid}"
    # Fallback if term_id missing
    term = str(df_row.get("term", "")).strip()
    return f"{cat}:{term}"


def top_terms_from_string(
    df: pd.DataFrame,
    *,
    top_n: int = 20,
    q_cut: float = 0.05,
    allowed_categories: Optional[Set[str]] = None,
    disallowed_categories: Optional[Set[str]] = None,
) -> List[str]:
    """
    Rank STRING enrichment by FDR (ascending), then p-value; return top term keys ≤ q_cut.
    Optionally filter by allowed/disallowed categories.
    """
    if df is None or df.empty:
        return []
    # Optional category filtering (e.g., keep only GO/KEGG/Reactome/WikiPathways)
    if "category" in df.columns:
        if allowed_categories:
            df = df[df["category"].isin(allowed_categories)].copy()
        if disallowed_categories:
            df = df[~df["category"].isin(disallowed_categories)].copy()
        if df.empty:
            return []
    # Normalize column names possibly present: 'fdr', 'p_value'
    cols = {c.lower(): c for c in df.columns}
    fdr_col = cols.get("fdr")
    p_col = cols.get("p_value") or cols.get("pvalue") or cols.get("p")
    work = df.copy()
    if fdr_col is None:
        # If FDR not provided, approximate by sorting p and treating p as q for cutoff
        work["_fdr"] = work[p_col] if p_col in work.columns else 1.0
        fdr_col = "_fdr"
    work = work.sort_values([fdr_col] + ([p_col] if p_col in work.columns else []), ascending=True)
    # filter by q_cut
    work = work[work[fdr_col] <= q_cut] if fdr_col in work.columns else work
    # build stable keys
    keys = [
        _string_term_key(row)
        for _, row in work.head(top_n).iterrows()
    ]
    return keys


def string_enrichment_ses_from_ensps(
    *,
    query_ensps_a: Sequence[str],
    query_ensps_b: Sequence[str],
    species: int = 9606,
    caller: str = "local_script",
    ensure_string_ids: bool = True,
    top_n: int = 20,
    q_cut: float = 0.05,
    return_tables: bool = True,
    allowed_categories: Optional[Set[str]] = None,
    disallowed_categories: Optional[Set[str]] = None,
    session: Optional[requests.Session] = None,
) -> Dict[str, Union[float, List[str], pd.DataFrame]]:
    """
    1) STRING API enrichment for A and B (starting from ENSP or STRING IDs)
    2) Rank by FDR and compute SES on top-N terms (≤ q_cut)
    Returns {'SES': float, 'ranked_terms_a': [...], 'ranked_terms_b': [...], 'enrich_a': df, 'enrich_b': df}
    """
    df_a = string_enrichment_from_ids(query_ensps_a, species=species, ensure_string_ids=ensure_string_ids, caller=caller, output_format="tsv", session=session)
    df_b = string_enrichment_from_ids(query_ensps_b, species=species, ensure_string_ids=ensure_string_ids, caller=caller, output_format="tsv", session=session)

    # category filters
    default_disallowed = {"PMID", "TISSUES"}
    use_allowed = allowed_categories
    use_disallowed = disallowed_categories if disallowed_categories is not None else default_disallowed

    ranked_a = top_terms_from_string(
        df_a, top_n=top_n, q_cut=q_cut,
        allowed_categories=use_allowed,
        disallowed_categories=use_disallowed,
    )
    ranked_b = top_terms_from_string(
        df_b, top_n=top_n, q_cut=q_cut,
        allowed_categories=use_allowed,
        disallowed_categories=use_disallowed,
    )

    N = min(len(ranked_a), len(ranked_b))
    ses = ses_from_ranked_terms(ranked_a, ranked_b, N=N) if N > 0 else 0.0

    out: Dict[str, Union[float, List[str], pd.DataFrame]] = {
        "SES": float(ses),
        "ranked_terms_a": ranked_a,
        "ranked_terms_b": ranked_b,
    }
    if return_tables:
        out.update({"enrich_a": df_a, "enrich_b": df_b})
    return out


###########################################
# 7) Per-category SES (GO/KEGG/Reactome/WikiPathways)
###########################################

def string_enrichment_ses_by_category(
    *,
    query_ensps_a: Sequence[str],
    query_ensps_b: Sequence[str],
    categories: Sequence[str],
    species: int = 9606,
    caller: str = "local_script",
    ensure_string_ids: bool = True,
    top_n: int = 20,
    q_cut: float = 0.05,
    return_tables: bool = True,
    disallowed_categories: Optional[Set[str]] = None,
    session: Optional[requests.Session] = None,
) -> Dict[str, Dict[str, Union[float, List[str], pd.DataFrame]]]:
    """
    Compute STRING enrichment once per set, then produce SES **per category** and an aggregated combined summary.
    Returns a dict:
      {
        'combined': {'SES': float, 'N_total': int, 'enrich_a': df, 'enrich_b': df},
        'by_category': {
            'Process': {'SES': ..., 'N': int, 'ranked_terms_a': [...], 'ranked_terms_b': [...]},
            'KEGG': {...}, ...
        }
      }
    Notes:
    - Combined SES is **not** computed by mixing terms across categories; it is the **weighted average** of per-category SES using weight N = min(|TopA|, |TopB|) per category.
    """
    # 1) Run enrichment once per side
    df_a = string_enrichment_from_ids(query_ensps_a, species=species, ensure_string_ids=ensure_string_ids, caller=caller, output_format="tsv", session=session)
    df_b = string_enrichment_from_ids(query_ensps_b, species=species, ensure_string_ids=ensure_string_ids, caller=caller, output_format="tsv", session=session)

    out: Dict[str, Dict[str, Union[float, List[str], pd.DataFrame]]] = {
        'combined': {},
        'by_category': {}
    }
    if return_tables:
        out['combined'].update({'enrich_a': df_a, 'enrich_b': df_b})

    # 2) Per-category
    for cat in categories:
        ra = top_terms_from_string(df_a, top_n=top_n, q_cut=q_cut, allowed_categories={cat})
        rb = top_terms_from_string(df_b, top_n=top_n, q_cut=q_cut, allowed_categories={cat})
        N = min(len(ra), len(rb))
        ses = ses_from_ranked_terms(ra, rb, N=N) if N > 0 else 0.0
        out['by_category'][cat] = {
            'SES': float(ses),
            'N': int(N),
            'ranked_terms_a': ra,
            'ranked_terms_b': rb,
        }

    # 3) Aggregated combined SES (weighted average, weights = N per category)
    Ns = [rec['N'] for rec in out['by_category'].values() if rec['N'] > 0]
    if Ns:
        weights = np.array(Ns, dtype=float)
        ses_vals = np.array([rec['SES'] for rec in out['by_category'].values() if rec['N'] > 0], dtype=float)
        agg = float(np.average(ses_vals, weights=weights))
        out['combined'].update({'SES': agg, 'N_total': int(np.sum(weights))})
    else:
        out['combined'].update({'SES': 0.0, 'N_total': 0})
    return out


###########################################
# 8) SES randomization test (empirical p)
###########################################

#
# ------------------------------
# Background from HUMAN PPI (SQLite)
# ------------------------------

def _to_ensp_core(x: str) -> str:
    """Normalize various STRING/Ensembl protein IDs to bare ENSP core.
    - '9606.ENSP00000354587.2' -> 'ENSP00000354587'
    - 'ENSP00000354587.2' -> 'ENSP00000354587'
    - Other prefixes are ignored.
    """
    if not isinstance(x, str) or not x:
        return ""
    s = x.strip()
    # strip STRING species prefix if present
    if s.startswith("9606."):
        s = s.split(".", 1)[1]
    # keep ENSP and drop version
    if s.startswith("ENSP"):
        return _normalize_ensembl_id(s)
    return ""


def _fetch_ppi_rows(
    db_path_local: str | Path | None = None,
    table: str = "STRING_PPI",
    col_u: str = "protein1",
    col_v: str = "protein2",
) -> List[Tuple[str, str]]:
    """Try to fetch (u,v) rows from a PPI table, with fallbacks for common table/column names.
    Returns a list of (u,v) strings; empty list if nothing found.
    """
    path = db_path_local or db_path
    try:
        con = sqlite3.connect(path)
        cur = con.cursor()
        # candidate tables & column-pairs to try
        table_candidates = [table, "Human_PPI", "STRING_PPI", "PPI", "ppi", "edges"]
        col_pairs = [(col_u, col_v), ("protein1","protein2"), ("node1","node2"), ("a","b"), ("source","target")]
        for t in table_candidates:
            # check table existence
            cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (t,))
            if not cur.fetchone():
                continue
            # try column pairs
            for u, v in col_pairs:
                try:
                    cur.execute(f"SELECT {u}, {v} FROM {t} LIMIT 1")
                except Exception:
                    continue
                # if above didn't error, fetch all rows
                cur.execute(f"SELECT {u}, {v} FROM {t}")
                rows = cur.fetchall()
                con.close()
                return rows or []
        con.close()
    except Exception:
        pass
    return []

def _load_background_from_ppi(
    db_path_local: str | Path | None = None,
    table: str = "STRING_PPI",
    col_u: str = "protein1",
    col_v: str = "protein2",
) -> List[str]:
    """Load unique ENSP IDs from a HUMAN PPI edge table (two protein columns).
    Assumes STRING-style IDs in columns (e.g., '9606.ENSP...') or ENSP.
    Returns a de-duplicated list of ENSP cores.
    """
    path = db_path_local or db_path
    try:
        rows = _fetch_ppi_rows(path, table=table, col_u=col_u, col_v=col_v)
        bg = set()
        for a, b in rows:
            ea = _to_ensp_core(a)
            eb = _to_ensp_core(b)
            if ea: bg.add(ea)
            if eb: bg.add(eb)
        return sorted(bg)
    except Exception:
        return []


def _ppi_degrees(
    db_path_local: str | Path | None = None,
    table: str = "STRING_PPI",
    col_u: str = "protein1",
    col_v: str = "protein2",
) -> Dict[str, int]:
    """Compute degrees for ENSP nodes from a PPI table."""
    path = db_path_local or db_path
    try:
        rows = _fetch_ppi_rows(path, table=table, col_u=col_u, col_v=col_v)
        deg: Dict[str, int] = {}
        for a, b in rows:
            ea = _to_ensp_core(a)
            eb = _to_ensp_core(b)
            if ea and eb:
                deg[ea] = deg.get(ea, 0) + 1
                deg[eb] = deg.get(eb, 0) + 1
        return deg
    except Exception:
        return {}


# === Inserted helper functions for degree binning and background loading ===
def _prepare_degree_bins(
    background_ids: Sequence[str],
    degrees: Dict[str, int],
    bins: int = 5,
) -> Tuple[np.ndarray, Dict[int, List[str]]]:
    """Precompute quantile bin edges and map background ids to bins.
    Returns (edges, bin2bg) where bin2bg[k] = list of background ids in bin k.
    """
    bg_deg = np.array([degrees.get(x, 0) for x in background_ids], dtype=float)
    if len(background_ids) < bins:
        bins = max(1, min(bins, len(background_ids)))
    qs = np.linspace(0, 1, bins + 1)
    edges = np.unique(np.quantile(bg_deg, qs))
    from collections import defaultdict
    bin2bg: Dict[int, List[str]] = defaultdict(list)
    if len(edges) > 1:
        def assign(vals):
            return np.digitize(vals, edges, right=True)
        bg_bins = assign(bg_deg)
        for x, b in zip(background_ids, bg_bins):
            bin2bg[int(b)].append(x)
    return edges, bin2bg


def _load_background_from_ppi_as_string_ids(
    db_path_local: str | Path | None = None,
    table: str = "STRING_PPI",
    col_u: str = "protein1",
    col_v: str = "protein2",
) -> List[str]:
    """Return unique STRING IDs (e.g., '9606.ENSP...') directly from a PPI table.
    Useful to skip get_string_ids during randomization (set ensure_string_ids=False for nulls).
    """
    rows = _fetch_ppi_rows(db_path_local or db_path, table=table, col_u=col_u, col_v=col_v)
    sids = set()
    for a, b in rows:
        if isinstance(a, str) and a.startswith("9606."):
            sids.add(a.strip())
        if isinstance(b, str) and b.startswith("9606."):
            sids.add(b.strip())
    return sorted(sids)


def _degree_matched_sample(
    target_ids: Sequence[str],
    background_ids: Sequence[str],
    degrees: Dict[str, int],
    bins: int = 5,
    rng: np.random.Generator | None = None,
) -> List[str]:
    """Sample |target| background nodes matching the degree histogram of target.
    - Assign background and target to quantile bins of degree.
    - For each bin, sample as many as the target has in that bin.
    Fallback: if a bin lacks enough candidates, sample from nearest bins.
    """
    if rng is None:
        rng = np.random.default_rng()
    # Build degree arrays
    bg_deg = np.array([degrees.get(x, 0) for x in background_ids], dtype=float)
    tgt_deg = np.array([degrees.get(_to_ensp_core(x), 0) for x in target_ids], dtype=float)
    # Define bin edges from background distribution (quantiles)
    if len(background_ids) < bins:
        bins = max(1, min(bins, len(background_ids)))
    qs = np.linspace(0, 1, bins + 1)
    edges = np.unique(np.quantile(bg_deg, qs))
    if len(edges) <= 1:
        # no variability; simple random sample
        return rng.choice(background_ids, size=len(target_ids), replace=False).tolist()
    # Assign to bins
    def assign(vals):
        return np.digitize(vals, edges, right=True)
    bg_bins = assign(bg_deg)
    tgt_bins = assign(tgt_deg)
    # Index background by bin
    from collections import defaultdict
    bin2bg = defaultdict(list)
    for x, b in zip(background_ids, bg_bins):
        bin2bg[int(b)].append(x)
    # Sample per bin
    out = []
    for b in np.unique(tgt_bins):
        need = int(np.sum(tgt_bins == b))
        pool = bin2bg.get(int(b), [])
        if len(pool) >= need:
            out.extend(rng.choice(pool, size=need, replace=False).tolist())
        else:
            # fallback: aggregate from neighboring bins
            take = []
            # try same bin first
            take.extend(pool)
            remain = need - len(take)
            # expand window
            step = 1
            while remain > 0 and (b - step >= 0 or b + step <= len(edges)):
                for nb in [b - step, b + step]:
                    if nb in bin2bg and nb >= 0:
                        cand = bin2bg[nb]
                        if cand:
                            k = min(remain, len(cand))
                            sel = rng.choice(cand, size=k, replace=False).tolist()
                            take.extend(sel)
                            remain -= k
                            if remain <= 0:
                                break
                step += 1
            # if still short, fill randomly from all background
            if remain > 0:
                rest = rng.choice(background_ids, size=remain, replace=False).tolist()
                take.extend(rest)
            out.extend(take[:need])
    # final size guard
    if len(out) > len(target_ids):
        out = rng.choice(out, size=len(target_ids), replace=False).tolist()
    return out

def _load_background_ensps_from_db(db_path_local: str | Path | None = None) -> List[str]:
    """
    Load a background of ENSP IDs from the local SQLite mapping table.
    Falls back to empty list if the table/column is missing.
    """
    path = db_path_local or db_path
    try:
        con = sqlite3.connect(path)
        cur = con.cursor()
        # Expect column 'ensp' with ENSP ids
        cur.execute("SELECT ensp FROM Ensembl_xref WHERE ensp IS NOT NULL")
        rows = cur.fetchall()
        con.close()
        bg = sorted({ _normalize_ensembl_id(r[0]) for r in rows if r and r[0] })
        return [x for x in bg if x]
    except Exception:
        return []


def string_enrichment_ses_randomized(
    *,
    query_ensps_a: Sequence[str],
    query_ensps_b: Sequence[str],
    n_iter: int = 1000,
    species: int = 9606,
    caller: str = "local_script",
    ensure_string_ids: bool = True,
    top_n: int = 20,
    q_cut: float = 0.05,
    allowed_categories: Optional[Set[str]] = None,
    disallowed_categories: Optional[Set[str]] = None,
    background_ensps: Optional[Sequence[str]] = None,
    db_path_local: Optional[str | Path] = None,
    background_source: str = "xref",   # 'xref' or 'ppi'
    ppi_table: str = "STRING_PPI",
    ppi_col_u: str = "protein1",
    ppi_col_v: str = "protein2",
    degree_matched: bool = False,
    degree_bins: int = 5,
    seed: Optional[int] = None,
    return_null: bool = False,
) -> Dict[str, Union[float, List[float]]]:
    """
    Empirical test of SES by randomizing both sets A and B with size-matched samples.

    Steps
    -----
    1) Compute observed SES using `string_enrichment_ses_from_ensps` with the provided filters.
    2) Build a background pool of ENSP IDs:
       - if `background_ensps` is given, use it;
       - else try to load from local SQLite `Ensembl_xref.ensp` or PPI.
       The background should correspond to the same species.
    3) Repeat n_iter times:
       - Sample |A| and |B| proteins without replacement from the background.
         (Degree-matched sampling if requested.)
       - Run STRING enrichment and SES with the same options.
       - Collect null SES values.
    4) Return observed SES, null mean/std, and empirical p = (1 + #null >= obs)/(n_iter + 1).
    """
    # 1) observed SES
    sess = requests.Session()
    try:
        observed = string_enrichment_ses_from_ensps(
            query_ensps_a=query_ensps_a,
            query_ensps_b=query_ensps_b,
            species=species,
            caller=caller,
            ensure_string_ids=ensure_string_ids,
            top_n=top_n,
            q_cut=q_cut,
            allowed_categories=allowed_categories,
            disallowed_categories=disallowed_categories,
            return_tables=False,
            session=sess,
        )
        ses_obs = float(observed.get("SES", 0.0))

        # 2) background pool
        if background_ensps is not None:
            bg = list(background_ensps)
        else:
            if background_source == "ppi":
                bg = _load_background_from_ppi(db_path_local=db_path_local, table=ppi_table, col_u=ppi_col_u, col_v=ppi_col_v)
            else:
                bg = _load_background_ensps_from_db(db_path_local)
        if not bg:
            raise ValueError("Background ENSP pool is empty. Provide `background_ensps` or a valid background source/table.")

        # Make sure inputs exist in background to avoid sampling bias
        bg = sorted({_normalize_ensembl_id(x) for x in bg if x})
        n_a = len({_normalize_ensembl_id(x) for x in query_ensps_a})
        n_b = len({_normalize_ensembl_id(x) for x in query_ensps_b})
        if n_a == 0 or n_b == 0:
            return {"SES": ses_obs, "p_emp": 1.0, "null_mean": float("nan"), "null_std": float("nan")}
        if len(bg) < max(n_a, n_b):
            raise ValueError(f"Background size ({len(bg)}) is smaller than required sample sizes (A={n_a}, B={n_b}).")

        # Pre-compute degree cache once if using PPI degree-matched sampling
        deg_cache: Dict[str, int] | None = None
        if degree_matched and background_source == "ppi":
            deg_cache = _ppi_degrees(db_path_local=db_path_local, table=ppi_table, col_u=ppi_col_u, col_v=ppi_col_v)

        rng = np.random.default_rng(seed)
        null_vals: List[float] = []

        for _ in range(int(max(1, n_iter))):
            if degree_matched and background_source == "ppi":
                samp_a = _degree_matched_sample(query_ensps_a, bg, deg_cache or {}, bins=degree_bins, rng=rng)
                samp_b = _degree_matched_sample(query_ensps_b, bg, deg_cache or {}, bins=degree_bins, rng=rng)
            else:
                samp_a = rng.choice(bg, size=n_a, replace=False).tolist()
                samp_b = rng.choice(bg, size=n_b, replace=False).tolist()
            try:
                res = string_enrichment_ses_from_ensps(
                    query_ensps_a=samp_a,
                    query_ensps_b=samp_b,
                    species=species,
                    caller=caller,
                    ensure_string_ids=ensure_string_ids,
                    top_n=top_n,
                    q_cut=q_cut,
                    allowed_categories=allowed_categories,
                    disallowed_categories=disallowed_categories,
                    return_tables=False,
                    session=sess,
                )
                null_vals.append(float(res.get("SES", 0.0)))
            except Exception:
                # Robust: skip failed iterations
                continue

        if not null_vals:
            return {"SES": ses_obs, "p_emp": 1.0, "null_mean": float("nan"), "null_std": float("nan")}

        null_arr = np.asarray(null_vals, dtype=float)
        null_mean = float(np.mean(null_arr))
        null_std  = float(np.std(null_arr, ddof=1)) if len(null_arr) > 1 else 0.0

        # empirical one-sided p-value (H0: SES_random >= SES_obs as strong or stronger)
        p_emp = (1.0 + float(np.sum(null_arr >= ses_obs))) / (len(null_arr) + 1.0)

        out: Dict[str, Union[float, List[float]]] = {
            "SES": ses_obs,
            "p_emp": float(p_emp),
            "null_mean": null_mean,
            "null_std": null_std,
        }
        if return_null:
            out["null_vals"] = null_vals
        return out
    finally:
        try:
            sess.close()
        except Exception:
            pass


###########################################
# 9) Randomization per-category + combined helper
###########################################

def string_enrichment_ses_randomized_by_category(
    *,
    query_ensps_a: Sequence[str],
    query_ensps_b: Sequence[str],
    combined_categories: Optional[Set[str]] = None,
    categories: Sequence[str] = ("Process","Function","Component","KEGG","Reactome","WikiPathways"),
    n_iter: int = 500,
    species: int = 9606,
    caller: str = "local_script",
    ensure_string_ids: bool = True,
    top_n: int = 20,
    q_cut: float = 0.05,
    db_path_local: Optional[str | Path] = None,
    background_source: str = "ppi",
    ppi_table: str = "STRING_PPI",
    ppi_col_u: str = "protein1",
    ppi_col_v: str = "protein2",
    degree_matched: bool = False,
    degree_bins: int = 5,
    seed: Optional[int] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Compute empirical SES (randomization) for:
      - a COMBINED set of categories (default: FUNCTIONAL_CATEGORIES), and
      - EACH category in `categories` individually.

    Returns mapping:
      {
        'Combined': { 'SES': float, 'p_emp': float, 'null_mean': float, 'null_std': float },
        '<cat>':   { 'SES': float, 'p_emp': float, 'null_mean': float, 'null_std': float }, ...
      }
    """
    if combined_categories is None:
        combined_categories = FUNCTIONAL_CATEGORIES

    sess = requests.Session()
    try:
        # -----------------------------
        # 0) Observed: compute once
        # -----------------------------
        df_a_obs = string_enrichment_from_ids(query_ensps_a, species=species, ensure_string_ids=ensure_string_ids, caller=caller, output_format="tsv", session=sess)
        df_b_obs = string_enrichment_from_ids(query_ensps_b, species=species, ensure_string_ids=ensure_string_ids, caller=caller, output_format="tsv", session=sess)

        def _percat_ses(df_a: pd.DataFrame, df_b: pd.DataFrame, cats: Sequence[str]) -> Tuple[Dict[str, float], Dict[str, int]]:
            ses_map: Dict[str, float] = {}
            n_map: Dict[str, int] = {}
            for cat in cats:
                ra = top_terms_from_string(df_a, top_n=top_n, q_cut=q_cut, allowed_categories={cat})
                rb = top_terms_from_string(df_b, top_n=top_n, q_cut=q_cut, allowed_categories={cat})
                N = min(len(ra), len(rb))
                ses_map[cat] = ses_from_ranked_terms(ra, rb, N=N) if N > 0 else 0.0
                n_map[cat] = N
            return ses_map, n_map

        # Observed per-category + combined (weighted by N)
        obs_ses_map, obs_n_map = _percat_ses(df_a_obs, df_b_obs, categories)
        Ns = np.array([obs_n_map[c] for c in categories if obs_n_map[c] > 0], dtype=float)
        if Ns.size > 0:
            ses_vals = np.array([obs_ses_map[c] for c in categories if obs_n_map[c] > 0], dtype=float)
            obs_combined_ses = float(np.average(ses_vals, weights=Ns))
        else:
            obs_combined_ses = 0.0

        # -----------------------------
        # 1) Null: reuse enrichment per iteration
        # -----------------------------
        bg_string_ids: Optional[List[str]] = None
        if background_source == "ppi":
            bg_core = _load_background_from_ppi(db_path_local=db_path_local, table=ppi_table, col_u=ppi_col_u, col_v=ppi_col_v)
            bg_string_ids = _load_background_from_ppi_as_string_ids(db_path_local=db_path_local, table=ppi_table, col_u=ppi_col_u, col_v=ppi_col_v)
        else:
            bg_core = _load_background_ensps_from_db(db_path_local)
        if not bg_string_ids and not bg_core:
            raise ValueError("Background ENSP pool is empty for randomized_by_category.")

        # Prepare degree cache/bins if requested
        deg_cache: Dict[str, int] | None = None
        edges_bins = None
        bin2bg = None
        if degree_matched and background_source == "ppi":
            deg_cache = _ppi_degrees(db_path_local=db_path_local, table=ppi_table, col_u=ppi_col_u, col_v=ppi_col_v)
            edges_bins, bin2bg = _prepare_degree_bins(bg_core, deg_cache or {}, degree_bins)

        # Sizes
        n_a = len({_normalize_ensembl_id(x) for x in query_ensps_a})
        n_b = len({_normalize_ensembl_id(x) for x in query_ensps_b})
        if n_a == 0 or n_b == 0:
            out = {"Combined": {"SES": obs_combined_ses, "p_emp": 1.0, "null_mean": float("nan"), "null_std": float("nan")}}
            for cat in categories:
                out[str(cat)] = {"SES": obs_ses_map.get(cat, 0.0), "p_emp": 1.0, "null_mean": float("nan"), "null_std": float("nan")}
            return out

        rng = np.random.default_rng(seed)
        null_per_cat: Dict[str, List[float]] = {str(c): [] for c in categories}
        null_combined: List[float] = []

        # Precompute core->STRING map if available
        core2sid: Dict[str, str] = {}
        if bg_string_ids:
            for sid in bg_string_ids:
                core = _to_ensp_core(sid)
                if core and core not in core2sid:
                    core2sid[core] = sid

        for _ in range(int(max(1, n_iter))):
            # Sample cores or STRING ids
            if background_source == "ppi":
                if degree_matched:
                    samp_a_core = _degree_matched_sample(query_ensps_a, bg_core, deg_cache or {}, bins=degree_bins, rng=rng, edges=edges_bins, bin2bg=bin2bg)
                    samp_b_core = _degree_matched_sample(query_ensps_b, bg_core, deg_cache or {}, bins=degree_bins, rng=rng, edges=edges_bins, bin2bg=bin2bg)
                else:
                    samp_a_core = rng.choice(bg_core, size=n_a, replace=False).tolist()
                    samp_b_core = rng.choice(bg_core, size=n_b, replace=False).tolist()
                # Map to STRING ids if available, else fall back to cores + ensure map
                if core2sid:
                    samp_a_ids = [core2sid.get(_normalize_ensembl_id(x), _normalize_ensembl_id(x)) for x in samp_a_core]
                    samp_b_ids = [core2sid.get(_normalize_ensembl_id(x), _normalize_ensembl_id(x)) for x in samp_b_core]
                    ensure_ids = False
                else:
                    samp_a_ids, samp_b_ids = samp_a_core, samp_b_core
                    ensure_ids = True
            else:
                # xref background (no STRING ids available)
                samp_a_ids = rng.choice(bg_core, size=n_a, replace=False).tolist()
                samp_b_ids = rng.choice(bg_core, size=n_b, replace=False).tolist()
                ensure_ids = True

            try:
                # One enrichment per side per iteration
                df_a = string_enrichment_from_ids(samp_a_ids, species=species, ensure_string_ids=ensure_ids, caller=caller, output_format="tsv", session=sess)
                df_b = string_enrichment_from_ids(samp_b_ids, species=species, ensure_string_ids=ensure_ids, caller=caller, output_format="tsv", session=sess)

                # Compute per-category SES from the same dfs
                ses_map, n_map = _percat_ses(df_a, df_b, categories)

                # Combined weighted by N
                Ns_it = np.array([n_map[c] for c in categories if n_map[c] > 0], dtype=float)
                if Ns_it.size > 0:
                    ses_vals_it = np.array([ses_map[c] for c in categories if n_map[c] > 0], dtype=float)
                    comb_it = float(np.average(ses_vals_it, weights=Ns_it))
                else:
                    comb_it = 0.0

                # Collect
                null_combined.append(comb_it)
                for c in categories:
                    null_per_cat[str(c)].append(float(ses_map.get(c, 0.0)))
            except Exception:
                continue

        # Summaries
        def _summ(obs: float, arr: List[float]) -> Dict[str, float]:
            if not arr:
                return {"SES": obs, "p_emp": 1.0, "null_mean": float("nan"), "null_std": float("nan")}
            a = np.asarray(arr, dtype=float)
            p_emp = (1.0 + float(np.sum(a >= obs))) / (len(a) + 1.0)
            return {"SES": obs, "p_emp": float(p_emp), "null_mean": float(np.mean(a)), "null_std": float(np.std(a, ddof=1)) if len(a) > 1 else 0.0}

        out: Dict[str, Dict[str, float]] = {}
        out["Combined"] = _summ(obs_combined_ses, null_combined)
        for c in categories:
            out[str(c)] = _summ(obs_ses_map.get(c, 0.0), null_per_cat[str(c)])
        return out
    finally:
        try:
            sess.close()
        except Exception:
            pass


def main():
    parser = argparse.ArgumentParser(description="STRING-based enrichment + SES (with optional per-category and randomization)")
    parser.add_argument("--herb-id", default="HERB002168", help="Herb ID to fetch ENSP list via Utils.get_ensp_ids")
    parser.add_argument("--disease-id", default="HBDIS001345", help="Disease ID to fetch ENSP list via Utils.get_ensp_ids")
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--q-cut", type=float, default=0.05)
    parser.add_argument("--categories", default="functional", choices=["go","pathways","functional","all"],
                        help="go=GO only; pathways=KEGG+Reactome+WP; functional=GO+KEGG+Reactome+WP; all=everything")
    parser.add_argument("--per-category", action="store_true", help="Report SES per category as well as combined")
    parser.add_argument("--randomize", type=int, default=0, help="#iterations for empirical p-value (0=skip)")
    parser.add_argument("--bg", default="ppi", choices=["xref","ppi"], help="Randomization background: xref table or PPI (default: Human PPI)")
    parser.add_argument("--degree-matched", action="store_true", help="Use degree-matched sampling when --bg=ppi")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ppi-table", default="Human_PPI", help="PPI table name (default: Human_PPI)")
    parser.add_argument("--ppi-col-u", default="protein1", help="PPI first protein column name")
    parser.add_argument("--ppi-col-v", default="protein2", help="PPI second protein column name")

    args = parser.parse_args()

    # Resolve categories preset
    if args.categories == "go":
        allowed = {"Process","Function","Component"}
    elif args.categories == "pathways":
        allowed = {"KEGG","Reactome","WikiPathways"}
    elif args.categories == "functional":
        allowed = FUNCTIONAL_CATEGORIES
    else:  # all
        allowed = None  # no filtering

    # Fetch ENSP lists
    herb_ensp = Utils.get_ensp_ids(args.herb_id)
    disease_ensp = Utils.get_ensp_ids(args.disease_id)

    # 1) Report SES
    if args.per_category:
        # Determine category set based on the preset
        if args.categories == 'go':
            cats = ["Process","Function","Component"]
        elif args.categories == 'pathways':
            cats = ["KEGG","Reactome","WikiPathways"]
        elif args.categories == 'functional':
            cats = sorted(FUNCTIONAL_CATEGORIES)
        else:  # 'all' -> use functional set to avoid mixing unrelated categories
            cats = sorted(FUNCTIONAL_CATEGORIES)

        cat_report = string_enrichment_ses_by_category(
            query_ensps_a=herb_ensp,
            query_ensps_b=disease_ensp,
            categories=cats,
            species=args.species,
            caller="proximity_pipeline",
            ensure_string_ids=True,
            top_n=args.top_n,
            q_cut=args.q_cut,
            return_tables=False,
        )
        print("[STRING][Combined] SES:", cat_report['combined'].get('SES'), "N_total=", cat_report['combined'].get('N_total'))
        for cat, rec in cat_report['by_category'].items():
            print(f"[STRING][{cat}] SES:", rec['SES'], "N=", rec['N'])
    else:
        # Non per-category: restrict to the selected preset to avoid cross-category mixing
        res = string_enrichment_ses_from_ensps(
            query_ensps_a=herb_ensp,
            query_ensps_b=disease_ensp,
            species=args.species,
            caller="proximity_pipeline",
            ensure_string_ids=True,
            top_n=args.top_n,
            q_cut=args.q_cut,
            allowed_categories=allowed,
        )
        print("[STRING] Top A:", res.get("ranked_terms_a", [])[:10])
        print("[STRING] Top B:", res.get("ranked_terms_b", [])[:10])
        print("[STRING] SES:", res.get("SES"))

    # 3) Optional randomization
    if args.randomize > 0:
        if args.per_category:
            rnd_cat = string_enrichment_ses_randomized_by_category(
                query_ensps_a=herb_ensp,
                query_ensps_b=disease_ensp,
                combined_categories=FUNCTIONAL_CATEGORIES,
                categories=( ["Process","Function","Component"] if args.categories == 'go' else
                            ["KEGG","Reactome","WikiPathways"] if args.categories == 'pathways' else
                            sorted(FUNCTIONAL_CATEGORIES) ),
                n_iter=args.randomize,
                species=args.species,
                caller="proximity_pipeline",
                ensure_string_ids=True,
                top_n=args.top_n,
                q_cut=args.q_cut,
                db_path_local=db_path,
                background_source=args.bg,
                ppi_table=args.ppi_table,
                ppi_col_u=args.ppi_col_u,
                ppi_col_v=args.ppi_col_v,
                degree_matched=args.degree_matched,
                seed=args.seed,
            )
            print("[STRING][Randomization][Combined] SES=", rnd_cat["Combined"]["SES"],
                  " p_emp=", rnd_cat["Combined"].get("p_emp", float('nan')),
                  " null_mean=", rnd_cat["Combined"].get("null_mean", float('nan')),
                  " null_std=", rnd_cat["Combined"].get("null_std", float('nan')))
            for cat in ( ["Process","Function","Component"] if args.categories == 'go' else
                         ["KEGG","Reactome","WikiPathways"] if args.categories == 'pathways' else
                         sorted(FUNCTIONAL_CATEGORIES) ):
                rec = rnd_cat[str(cat)]
                print(f"[STRING][Randomization][{cat}] SES= {rec['SES']}  p_emp= {rec['p_emp']}  null_mean= {rec['null_mean']}  null_std= {rec['null_std']}")
        else:
            rnd = string_enrichment_ses_randomized(
                query_ensps_a=herb_ensp,
                query_ensps_b=disease_ensp,
                n_iter=args.randomize,
                species=args.species,
                caller="proximity_pipeline",
                ensure_string_ids=True,
                top_n=args.top_n,
                q_cut=args.q_cut,
                allowed_categories=allowed,
                db_path_local=db_path,
                background_source=args.bg,
                ppi_table=args.ppi_table,
                ppi_col_u=args.ppi_col_u,
                ppi_col_v=args.ppi_col_v,
                degree_matched=args.degree_matched,
                seed=args.seed,
                return_null=False,
            )
            print("[STRING][Randomization] SES=", rnd["SES"], " p_emp=", rnd["p_emp"], " null_mean=", rnd["null_mean"], " null_std=", rnd["null_std"])

if __name__ == "__main__":
    main()

