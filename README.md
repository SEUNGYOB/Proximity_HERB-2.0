# 🧬 Proximity-Final

Network-based **Herbal Medicine – Target – Disease** relationship analysis using **STRING Human PPI** and **HERB 2.0 (HIT)** database.  
Implements efficient **Network Proximity** (A–B set distance & Z-score) computation with **incremental memmap caching**.

---

## 📁 1. Database Download

```bash
python DB_download.py
```

[📥 Download Data (Google Drive)](https://drive.google.com/file/d/1MceeVsg9utgvFRjC73QPgJWPh9o9J_US/view?usp=drive_link)

This command:
- Creates the `Data/` folder (if missing)  
- Downloads the SQLite database (`DB.db`) from Google Drive  
- Verifies integrity and saves it locally

---

## 🧩 2. Database Source

### 🔹 2.1 Human PPI (STRING)

- **Source:** [STRING v12.0](https://string-db.org/cgi/download?sessionId=bxgqh7YBMvBD&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt)  
- **File:** `9606.protein.links.v12.0.txt.gz`  
- **Processed with:** `build_ppi_db()`  
- **Uploaded Table:** `Human_PPI`  
- **Total edges:** 473,860  
- **Node format:** `9606.ENSPxxxx` → prefix `9606.` removed during upload  

---

### 🔹 2.2 HERB 2.0 (HIT)

- **Source:** [HERB 2.0 Download Portal](http://47.92.70.12/Download/)  
- **Structure:** Formula–Herb–Ingredient–Target–Disease network  
- **Filtering:** ADME-passing ingredients only  
- **Target ID format:** `ENSG` (converted to `ENSP` via BioMart)

#### ⚗️ ADME Screening Criteria (Balanced Profile)

| Property | Threshold | Description |
|-----------|------------|--------------|
| OB_score | ≥ 30 | Oral bioavailability |
| Drug_likeness | ≥ 0.18 | QED-like metric |
| MolWt | ≤ 500 | Molecular weight |
| NumHAcceptors | ≤ 10 | H-bond acceptors |
| NumHDonors | ≤ 5 | H-bond donors |
| MolLogP | ≤ 5.5 | Lipophilicity |
| NumRotatableBonds | ≤ 10 | Flexibility |

---

## 📊 3. Database Overview (`DB.db`)

SQLite database containing curated HERB–STRING integration.  
Each table’s purpose and identifiers are summarized below.

| Table | Source | Description | Key Columns |
|--------|--------|-------------|--------------|
| **Ensembl_xref** | BioMart-derived | Canonical ENSG–ENSP mapping | `ensg`, `ensp`, `hgnc_symbol`, `uniprot` |
| **Human_PPI** | STRING v12.0 | Protein–protein interaction edges | `protein1`, `protein2` (ENSP) |
| **Herb_info** | HERB 2.0 | Herb metadata | `Herb_id`, names, categories |
| **Ingredient_info** | HERB 2.0 | Ingredient + ADME data | `Ingredient_id`, `Canonical_smiles`, ADME fields |
| **Target_info** | HERB 2.0 | Target genes/proteins | `Target_id`, `Ensembl_id` |
| **Disease_info** | HERB 2.0 | Disease metadata | `Disease_id`, `Disease_name` |
| **HERB_edges** | HERB 2.0 | Full heterogeneous edge list | `Node_1`, `Node_2` |

> **Workflow Summary:**  
> - Map `Target (ENSG)` → `Protein (ENSP)` using `Ensembl_xref`  
> - Integrate into `Human_PPI` for network proximity computation  
> - PPI nodes are in `ENSP` format (`9606.` prefix removed)

---

## ⚙️ 4. Command-Line Usage — `run_proximity.py`

CLI for computing **network proximity (Z-score)** between two biological node sets  
(e.g., Herb–Disease, Ingredient–Target) using the HERB–STRING Human PPI network.

### 🔧 Basic Usage

```bash
python run_proximity.py --a-id HERB002168 --b-id HBDIS001345
```

### 🧠 Optional Arguments

| Argument | Default | Description |
|-----------|----------|-------------|
| `--db` | `./Data/DB.db` | Path to SQLite database |
| `--table` | `Human_PPI` | Edge table name |
| `--src-col` | `protein1` | Source column name |
| `--tgt-col` | `protein2` | Target column name |
| `--a-id` | — | Source node-set ID |
| `--b-id` | — | Target node-set ID |
| `--random-time` | 100 | # of degree-matched random samples |
| `--seed` | 42 | Random seed |
| `--workers` | 16 | Number of worker threads |
| `--out` | — | Optional output JSON file |

### 💡 Example Commands

```bash
# Basic
python run_proximity.py --a-id HERB002168 --b-id HBDIS001345

# More iterations
python run_proximity.py --a-id HERB002168 --b-id HBDIS001345 --random-time 200 --workers 32

# Save results
python run_proximity.py --a-id HERB002168 --b-id HBDIS001345 --out ./Results/HERB002168__HBDIS001345.json
```

### 📊 Example Output

```json
{
  "A_id": "HERB002168",
  "B_id": "HBDIS001345",
  "random_time": 100,
  "seed": 42,
  "result": {
    "shortest": 1.7488,
    "Z_score": { "d": 1.7488, "z": -7.79, "mean": 2.0093, "std": 0.0334, "p": 0.0 }
  }
}
```

> **Note:**  
> - Distance caching uses incremental memmap (`./Data/Human_PPI/rows_mm.npy`)  
> - Lightweight mode — no full adjacency matrix stored  
> - Suitable for batch runs and Snakemake workflows

---

# 🧩 Evaluation.py — Functional Proximity & Enrichment Analysis

`Evaluation.py` assesses **functional similarity between Herb and Disease target sets**  
based on STRING enrichment and **Shared Enrichment Score (SES)**.

---

## 🚀 Main Features
- STRING Enrichment: Fetches GO/Pathway enrichment automatically  
- Shared Enrichment Score (SES): Measures overlap in top enriched terms  
- Randomization Test: Computes empirical p-values (supports degree matching)  
- Per-Category Analysis: Separate SES for GO, KEGG, Reactome, WikiPathways  

---

## ⚙️ CLI Example

```bash
python Evaluation.py \
  --herb-id HERB002168 \
  --disease-id HBDIS001345 \
  --categories functional \
  --per-category \
  --randomize 10 \
  --seed 42
```

---

## 🧠 Key Options

| Option | Description | Default |
|--------|--------------|----------|
| `--herb-id` | HERB Herb ID | `HERB002168` |
| `--disease-id` | HERB Disease ID | `HBDIS001345` |
| `--categories` | Analysis scope (`go`, `pathways`, `functional`, `all`) | `functional` |
| `--per-category` | Report SES per category | `False` |
| `--randomize` | Randomization iterations (0 = skip) | `0` |
| `--bg` | Background source (`xref`, `ppi`) | `ppi` |
| `--degree-matched` | Enable degree-matched sampling | `False` |
| `--ppi-table` | PPI table name | `Human_PPI` |
| `--seed` | Random seed | `42` |

---

## 📈 Example Output

```
[STRING][Combined] SES: 0.1167  N_total=80  
[STRING][Function] SES: 0.0294  N=20  
[STRING][KEGG] SES: 0.1446  N=20  
[STRING][Process] SES: 0.1371  N=20  
[STRING][WikiPathways] SES: 0.1555  N=20  

[STRING][Randomization][Combined] SES= 0.1166  p_emp= 0.09  null_mean= 0.0  null_std= 0.0
```

---

## 🧬 Analysis Workflow
1. Load target protein ENSPs for Herb and Disease  
2. Run STRING enrichment (GO/KEGG/Reactome/WikiPathways)  
3. Compute Shared Enrichment Score (SES) from top overlapping terms  
4. (Optional) Run Randomization Test to estimate empirical significance (`p_emp`)  
5. (Optional) Report category-wise SES and combined weighted mean  

---

## 📦 Dependencies
- Python 3.9+  
- `requests`, `pandas`, `numpy`, `sqlite3`  
- Requires local SQLite: `Data/DB.db`  
- Internet access for STRING API

---

## 🔖 Notes
- Combined SES = weighted average (by N) of per-category SES  
- Default randomization background = Human PPI (`--bg ppi`)  
- Degree-matched sampling only applies to PPI-based background
