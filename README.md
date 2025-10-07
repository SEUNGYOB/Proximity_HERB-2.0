# 🧬 Proximity-Final

Network-based **Herbal Medicine – Target – Disease** relationship analysis using **STRING Human PPI** and **HERB 2.0 (HIT)** database.  
Implements efficient **Network Proximity** (A–B set distance & Z-score) calculation with **incremental memmap caching**.

---

## 📁 1. Database Download

### 🔹 1.1 Human PPI (STRING)
- **Source:** [STRING v12.0](https://string-db.org/cgi/download?sessionId=bxgqh7YBMvBD&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt)
- **File:** `9606.protein.links.v12.0.txt.gz`
- **Processed with:** `build_ppi_db()`
- **Uploaded Table:** `Human_PPI`
  - Total edges: **473,860**
  - Node format: `9606.ENSPXXXX` → prefix `9606.` removed during DB upload

---

### 🔹 1.2 HERB 2.0 (HIT) Data
- **Source:** [HERB 2.0 Download Portal](http://47.92.70.12/Download/)
- **Categories:**
  | Type | Range | Example ID | Note |
  |------|--------|-------------|------|
  | Herb | `HERB000001` ~ `HERB007263` | — |  |
  | Ingredient | `HBIN000001` ~ `HBIN049258` | — | ADME filtering applied |
  | Formula | `HBFO000001` ~ `HBFO006743` | — |  |
  | Target | `HBTAR000001` ~ `HBTAR015515` | ENSG format |  |
  | Disease | `HBDIS000001` ~ `HBDIS030170` | — |  |
  | Meta-analysis | `HBMA000001` | — |  |
  | Clinical trial | `HBCT000001` | — |  |
  | Reference | `HBREF000001` | — |  |
  | Experiment | `HBEXP000002` | — |  |

---

## 🧩 2. Database Construction

**Tables**
- `Human_PPI`
- `Herb_2_0_edge`
- `Herb_2_0_node`

**Filtering**
- Retained edges:  
  `Formula–Herb–Ingredient–Target (F–H–I–T)` and `Disease–Target (D–T)`
- Only **ADME-passing ingredients** are kept.

### ⚗️ ADME Screening Criteria (Balanced Profile)

| Property | Threshold | Rule |
|-----------|------------|------|
| OB_score | ≥ 30 | None = neutral |
| Drug_likeness | ≥ 0.18 |  |
| MolWt | ≤ 500 |  |
| NumHAcceptors | ≤ 10 |  |
| NumHDonors | ≤ 5 |  |
| MolLogP | ≤ 5.5 |  |
| NumRotatableBonds | ≤ 10 |  |

---
## 📊 Database Information (DB.db)`

SQLite `DB.db`의 주요 테이블을 워크북 시트로 내보낸 파일입니다.  
각 시트의 출처·의미·주요 식별자를 아래에 정리했습니다. *(괄호 안은 현재 행 수 예시)*

| Table | Source | What it contains | Key IDs / Notes |
|---|---|---|---|
| **BIOMART_mapping_raw** (533,785) | Ensembl **BioMart** | BioMart에서 그대로 받은 원시 매핑 | `Gene stable ID(ENSG)`, `Protein stable ID(ENSP)`, `Ensembl Canonical`, `APPRIS`, `Gene name` |
| **Ensembl_xref** (23,869) | Derived (from BioMart) | **대표 1:1 매핑** (ENSG→ENSP, Canonical/APPRIS 우선 규칙 적용) | `ensg`, `ensp`, `hgnc_symbol`, `uniprot`, `picked_by` |
| **Ensembl_xref_all** (245,535) | Derived (from BioMart) | **전체 1:N 매핑** (모든 ENSG–ENSP 쌍) | `ensg`, `ensp`, `hgnc_symbol`, `uniprot` |
| **Human_PPI** (473,860) | **STRING v12.0** | 인간 단백질–단백질 상호작용 엣지 | `protein1`, `protein2` (형태는 `ENSP`; 업로드 시 `9606.` 접두사 제거) |
| **Disease_info** (30,170) | **HERB 2.0 (HIT)** | 질병 메타데이터 | `Disease_id (HBDIS...)`, `Disease_name`, `DisGeNET_id`, 분류 체계/식별자 |
| **Formula_info** (6,743) | HERB 2.0 | 처방(방제) 메타데이터 | `Formula_id (HBFO...)`, 한/영/병음명, 범주/출전 |
| **Herb_info** (6,892) | HERB 2.0 | 본초(약재) 메타데이터 | `Herb_id (HERB...)`, 한/영/병음/이명 |
| **Ingredient_info** (44,595) | HERB 2.0 | 성분(화합물) 메타데이터 + **ADME 필드** | `Ingredient_id (HBIN...)`, `Canonical_smiles`, `Molecular_formula`, `OB_score`, `Drug_likeness`, `MolWt`, `NumHAcceptors/Donors`, `MolLogP`, `NumRotatableBonds` |
| **Target_info** (15,515) | HERB 2.0 | 타깃 메타데이터 | `Target_id (HBTAR...)`, `Ensembl_id (ENSG)` |
| **HERB_edges** (859,757) | **Scraped from HERB 2.0** | HERB 2.0에서 수집한 **원시 네트워크 엣지** | `Node_1`, `Node_2` (노드에는 `HERB/HBIN/HBTAR/HBDIS/HBFO` 등 코드가 혼재) |

> 🧭 **파이프라인 요약**  
> - **BioMart**로 ENSG→ENSP 매핑을 만들고(`Ensembl_xref`, `Ensembl_xref_all`),  
> - **HERB_edges**의 Target(주로 ENSG)을 **ENSP**로 연결해 **STRING Human_PPI** 네트워크에서 **Network Proximity**를 계산합니다.  
> - PPI 노드 표기는 `ENSP`; STRING 원본의 `9606.ENSPxxxxx`는 업로드 시 `9606.` 접두사를 제거했습니다.

### Export 재현
아래 스니펫으로 DB를 엑셀로 다시 덤프할 수 있습니다.
```python
# notebook_empty.py 내 함수 사용
export_db_to_excel("./Data/DB.db", "./HERB_KG_dump.xlsx", include_views=False)


## 🧬 3. BioMart (ENSG → ENSP Mapping)

- **Portal:** [Ensembl BioMart](https://www.ensembl.org/biomart/martview/8403dac70986fcef75d758c4bc648d6f)
- **Dataset:** *Ensembl Genes 115, Human (GRCh38.p14)*
- **Columns Used:**
  - Gene stable ID (ENSG)
  - Protein stable ID (ENSP)
  - APPRIS annotation
  - Ensembl Canonical
  - Gene name
- **Output:** `ENSG → ENSP` mapping table  
  Used to link HERB Target (ENSG) to STRING PPI (ENSP).

---

## 🔗 4. Network Proximity Calculation

- **Graph:** Human_PPI (STRING v12.0) — ENSP–ENSP edges  
- **Distance options:**
  - Precomputed full matrix (`DistMatrix`)  
  - Incremental memmap row caching (`rows_mm.npy`)

### 🧠 Output Format

```python
{
  "shortest": <observed mean distance>,
  "Z_score": {
      "d": <observed>,
      "z": <Z-score>,
      "mean": <random mean>,
      "std": <random std>,
      "p": <empirical p-value>
  }
}
