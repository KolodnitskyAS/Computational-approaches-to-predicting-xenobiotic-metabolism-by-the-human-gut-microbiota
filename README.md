#  Merged MASI–MDIPID Dataset

This folder provides the complete data-integration pipeline, analysis code, and documentation for the merged database analysis of microbiota-mediated xenobiotic metabolism records discussed in the manuscript. All scripts are provided in Python (3.x) and require pandas, numpy, and openpyxl.

## Contents

```
├── data
│   ├── MASI_MMDR_combined.xlsx              — combined dataset before enrichment (Step 3 output, 10,969 rows)
│   ├── MASI_MMDR_combined_enriched.xlsx     — final enriched dataset (pipeline output, 8,886 rows × 38 columns)
│   └── Enriched_Dataset.xlsx                — final dataset used for statistical analysis (8,886 rows × 38 columns)
├── full_pipeline.py                         — end-to-end 12-step merge pipeline (Steps 3–12)
├── statistical_analysis.py                  — reproduces all quantitative results reported in the manuscript
├── data_dictionary.csv                      — column-by-column description of the final dataset
├── statistical_analysis_summary.csv         — output of statistical_analysis.py
└── README.md                                — this file
```

## Database versions and access dates

MASI version 20200928
MDIPID (MMDR) current at access, 15 March 2025



## Pipeline overview (12 steps)

The full pipeline is implemented in full_pipeline.py and transforms two independently curated databases into a single, deduplicated, and enriched dataset. Steps 1–2 produce the pre-merged intermediate files (MASI_merged_enriched.xlsx and MMDR_merged_enriched.xlsx), which are loaded as input by the pipeline script. Steps 3–12 are executed sequentially within full_pipeline.py. The data flow is summarized below:

```
MASI_merged_enriched.xlsx (4,295 rows) ──┐
                                         ├──► MASI_MMDR_combined.xlsx (10,969 rows)
MMDR_merged_enriched.xlsx (6,674 rows) ──┘              │
                                                         ▼
                                          Internal enrichment (10,669 rows)
                                                         │
                                                         ▼
                                          Duplicate row merge (9,467 rows)
                                                         │
                                                         ▼
                                          Entity deduplication (9,420 rows)
                                                         │
                                                         ▼
                                          Small-molecule filter (8,886 rows)
                                                         │
                                                         ▼
                                          MASI_MMDR_combined_enriched.xlsx
```

### Step 1: Merge MASI main + supplementary files

The MASI main interaction table is joined with its supplementary taxonomy file (genus, family, taxonomic level) and probiotic annotation file (probiotic status, use case, research stage) via left joins on Microbe-Tax-ID / microbe_id. Column names are standardized (e.g., Microbe-Name → Microbe_Name, Substance-Name → Substance_Name, id_pubchem → PubChem_CID). Output: MASI_merged_enriched.xlsx (4,295 rows × 47 columns).

### Step 2: Merge MMDR main + supplementary files

The MMDR main interaction table is joined with three supplementary files: drug properties (PubChem CID, InChI, SMILES, physicochemical descriptors), protein annotations (EC classification, KEGG pathway, PDB ID, FASTA sequence, function description), and taxonomy (phylum, class, order, family, kingdom). Joins are performed on Drugid, Related_protein_id, and SpeciesID, respectively. Output: MMDR_merged_enriched.xlsx (6,674 rows × 91 columns).

### Step 3: Combine MASI + MMDR with cross-enrichment

The two enriched dataframes are combined using an outer-join strategy. A Source_Dataset column is added (values: MASI, MMDR, or MASI+MMDR). Overlapping microbe–substance pairs present in both databases (285 pairs identified) are cross-enriched: where one database has a missing field and the other provides a value, the missing cell is filled. Overlapping MMDR rows are absorbed into the corresponding MASI rows to avoid duplication. Output: MASI_MMDR_combined.xlsx (10,969 rows × 119 columns; 1,013 fields filled by cross-enrichment).

### Step 4: Internal enrichment (three levels)

Missing values within the combined dataset are propagated at three hierarchical levels:

- Substance level: If the same substance (by name) appears in multiple rows, non-ambiguous identifier values (PubChem CID, InChIKey, DrugBank ID, molecular formula, etc.) are propagated to rows where they are missing.
- Microbe level: Taxonomic and probiotic fields are propagated across rows sharing the same microbe name.
- Pair level (microbe + substance): Interaction-specific fields (metabolism type, mechanism, outcome, protein, EC number, kinetic parameters) are propagated where a pair has multiple rows and only one unique non-null value exists for a given field.

Only cells with a single unambiguous source value are filled; ambiguous cases (multiple conflicting values) are left as NaN. Result: 105,849 additional fields filled; 10,669 rows × 121 columns.

### Step 5: Audit and revert incorrectly filled pair-level cells

A quality-control step compares the enriched dataset against the pre-enrichment version (MASI_MMDR_combined.xlsx) to detect cells that were incorrectly filled at the pair level. Specifically, if a microbe–substance pair has multiple rows describing genuinely different enzymatic reactions (different proteins), protein-specific fields (Related_protein, EC number, UniProt, KEGG, FASTA, etc.) must not be copied across rows. Any such erroneously propagated values are reverted to NaN. Result: 197 cells reverted.

### Step 6: Flag substance and taxonomy conflicts

Two types of data conflicts are detected and flagged:

- Substance conflicts (31 substances): The same substance name maps to different chemical identifiers (PubChem CID, InChIKey, SMILES, ChEBI) across MASI and MMDR, indicating possible differences in salt form, stereoisomer, or annotation error.
- Taxonomy conflicts (12 microbes): The same microbe name is assigned different genus or family classifications in the two databases, reflecting outdated vs. current nomenclature.

Flag columns (Substance_Data_Conflict, Microbe_Taxonomy_Conflict) are added for traceability.

### Step 7: Unify taxonomy to modern nomenclature

Taxonomy conflicts identified in Step 6 are resolved by adopting current nomenclature (MMDR is more up-to-date). 

Result: 12 microbes updated, 85 cells changed.

### Step 8: Merge duplicate and semantically equivalent columns

Redundant or semantically overlapping columns are consolidated to reduce dimensionality:

- GeneName_Standarized + Related_protein_gene_name → Gene_Info (format: "gene_name | protein_description")
- ECID_Stand + EC1Name + EC2Name + EC3Name → EC_Classification (format: "2.7.1.95 | Transferase → Kinase → Phosphotransferase")
- Metabolites + Metabolites_CID → Metabolites_Info (format: "Metabolite_name (CID:12345)")

Conflict flag columns from Step 6 are removed. Result: 44 columns → 36 columns.

### Step 9: Merge duplicate interaction rows

Rows describing the same interaction (same microbe + substance + protein) are merged into a single record. Grouping logic:

- If all rows for a microbe–substance pair share the same protein annotation (or all lack one), they are merged into one row.
- If a pair has rows with genuinely different proteins (different enzymatic reactions), those rows are kept separate.
- Rows lacking protein annotation are absorbed into the single-protein group when only one protein exists for that pair.

Merge strategy for text fields: unique values are combined with " | " separator; reference IDs are concatenated with "; ". Result: 10,669 → 9,467 rows (1,202 rows eliminated); 120 pairs retain multiple rows representing genuinely distinct interactions.

### Step 10: Entity deduplication

Two types of entity-level deduplication are performed:

- Substance deduplication: Different substance names sharing the same PubChem CID are renamed to the most frequently occurring canonical name. After renaming, any newly created duplicate rows are re-merged. Result: 65 substance rows renamed; 1,319 → 1,251 unique substances.
- Microbe deduplication: Microbe names differing only in whitespace or case are normalized. Result: 815 → 814 unique microbes.

Overall: 9,467 → 9,420 rows.

### Step 11: Restore kinetic parameters (Km, Km_combination)

The KM (Michaelis constant) and KM_combination columns are restored from the original MMDR source data (`MMDR_merged_enriched.xlsx`). Matching is performed by microbe name + substance name + protein identity (where available). MASI records do not contain kinetic data and retain NaN for these fields. Result: 116 KM values and 56 KM_combination values restored.

### Step 12: Filter to small molecules (valid PubChem CID required)

Records lacking a valid PubChem CID are removed, as the dataset focuses on small-molecule drug metabolism. Excluded records include immunotherapy agents (anti-PD-1, anti-CTLA-4), biologics (infliximab, adalimumab), radiotherapy, drug class labels (NSAIDs, PPIs), and dietary compounds without chemical identifiers. Result: 9,420 → 8,886 rows (533 non-small-molecule records removed). Output: MASI_MMDR_combined_enriched.xlsx.

---

## Identifier standardization summary

- Compounds: Standardized to PubChem Compound ID (CID). Records for which no single-compound CID could be resolved were excluded in Step 12. 100% of retained records have a resolved PubChem CID. Additional cross-references include DrugBank, ChEBI, KEGG, and TTD identifiers.
- Microorganisms: Names harmonized to current taxonomic nomenclature (Step 7), reconciling superseded and current genus/species names across the two source databases.
- Enzymes: Annotated with EC classification numbers, UniProt accession numbers, gene names, and KEGG pathway assignments where resolvable (Steps 2, 4, 8).

## Deduplication summary

Duplicate records were identified and resolved at multiple levels (Steps 9–10) using composite keys: microbe name + substance name + protein annotation. Records sharing the same compound and microbe but attributed to distinct enzymes represent distinct biochemical observations and were retained separately. The manuscript-reported figure of 8,886 records corresponds to the dataset after all deduplication and filtering steps, not the raw merged count.

## Kinetic parameter handling

Km, Vmax, kcat, and CLint values, where reported in the source literature, were retained as structured entries (columns KM, KM_combination) together with their originally reported units (predominantly µM, mM, nM, or enzyme-activity units such as nmol·min⁻¹·mg protein⁻¹). No unit conversion or imputation of missing values was performed. The statistical_analysis.py script quantifies kinetic annotation coverage and reports the severe deficit of model-ready kinetic data as a principal barrier to PBPK modeling.

