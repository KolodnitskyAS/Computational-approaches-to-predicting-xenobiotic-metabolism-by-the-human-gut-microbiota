
"""
==========================================================================
Statistical Analysis of the Merged MASI–MDIPID Drug–Microbe Interaction
Dataset




  1. Total number of drug–microbe interaction records
  2. Experimental system breakdown (in vitro, in vivo, etc.)
  3. Methodological bias — HTS dominance among annotated records
  4. Protein-level annotation coverage
  5. EC classification coverage
  6. Michaelis constant (Km) coverage
  7. Km + combination kinetic parameter coverage
  8. Number of unique annotated enzymes with kinetic data
  9. Number and fraction of unique drugs with kinetic data
 10. Source database distribution (MASI / MMDR / MASI+MMDR)
 11. Microbe diversity statistics
 12. Drug diversity statistics

Input:  Enriched_Dataset.xlsx (8,886 rows x 38 columns)
Output: Printed statistical report + CSV summary table

Usage:
    python statistical_analysis.py
==========================================================================
"""

import pandas as pd
import numpy as np
import os
import sys

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATA_DIR = "/home/z/my-project/download"
INPUT_FILE = "Enriched_Dataset.xlsx"
OUTPUT_DIR = "/home/z/my-project/download"


def load_dataset():
    """Load the final enriched dataset."""
    path = os.path.join(DATA_DIR, INPUT_FILE)
    if not os.path.exists(path):
        # Try upload dir as fallback
        path = os.path.join("/home/z/my-project/upload", INPUT_FILE)
    if not os.path.exists(path):
        print(f"ERROR: Cannot find {INPUT_FILE}")
        sys.exit(1)
    df = pd.read_excel(path)
    print(f"Loaded: {path}")
    print(f"Shape: {df.shape}")
    return df


# ---------------------------------------------------------------------------
# 1. Total records
# ---------------------------------------------------------------------------
def stat_total_records(df):
    """Total number of drug–microbe interaction records."""
    n = len(df)
    print(f"\n{'='*70}")
    print(f"1. TOTAL RECORDS")
    print(f"{'='*70}")
    print(f"   Total drug–microbe interaction records: {n:,}")
    return {"metric": "Total records", "value": n, "pct": None}


# ---------------------------------------------------------------------------
# 2. Experimental system breakdown
# ---------------------------------------------------------------------------
def stat_experiment_system(df):
    """Breakdown of experimental systems (in vitro, in vivo, etc.)."""
    print(f"\n{'='*70}")
    print(f"2. EXPERIMENTAL SYSTEM BREAKDOWN")
    print(f"{'='*70}")

    total = len(df)
    col = "Experiment_System"

    # Classify each record
    is_in_vitro = df[col].str.contains(r"(?i)in\s*vitro", na=False)
    is_in_vivo = df[col].str.contains(r"(?i)in\s*vivo", na=False)
    is_human = df[col].str.contains(r"(?i)human|fecal|serum|urine|blood", na=False)
    is_empty = df[col].isna() | (df[col].astype(str).str.strip().isin(["", ".", "n.a."]))

    n_in_vitro = is_in_vitro.sum()
    n_in_vivo = is_in_vivo.sum()
    n_human = is_human.sum()
    n_empty = is_empty.sum()
    n_annotated = total - n_empty

    # Pure in vitro (not also in vivo)
    n_pure_in_vitro = (is_in_vitro & ~is_in_vivo).sum()

    print(f"   In vitro (any):          {n_in_vitro:>6}  ({n_in_vitro/total*100:.1f}%)")
    print(f"   In vitro (pure):         {n_pure_in_vitro:>6}  ({n_pure_in_vitro/total*100:.1f}%)")
    print(f"   In vivo (any):           {n_in_vivo:>6}  ({n_in_vivo/total*100:.1f}%)")
    print(f"   Human/fecal samples:     {n_human:>6}  ({n_human/total*100:.1f}%)")
    print(f"   Not annotated:           {n_empty:>6}  ({n_empty/total*100:.1f}%)")
    print(f"   Annotated:               {n_annotated:>6}  ({n_annotated/total*100:.1f}%)")

    # Detailed value counts for all system types
    print(f"\n   Detailed Experiment_System value counts:")
    sys_counts = df[col].value_counts()
    for val, cnt in sys_counts.items():
        print(f"     {val}: {cnt} ({cnt/total*100:.1f}%)")

    return [
        {"metric": "In vitro records", "value": n_pure_in_vitro,
         "pct": f"{n_pure_in_vitro/total*100:.1f}%"},
        {"metric": "In vivo records", "value": n_in_vivo,
         "pct": f"{n_in_vivo/total*100:.1f}%"},
        {"metric": "System annotated records", "value": n_annotated,
         "pct": f"{n_annotated/total*100:.1f}%"},
    ]


# ---------------------------------------------------------------------------
# 3. Methodological bias — HTS dominance
# ---------------------------------------------------------------------------
def stat_hts_bias(df):
    """Quantify the dominance of high-throughput screening among
    methodologically annotated records."""
    print(f"\n{'='*70}")
    print(f"3. METHODOLOGICAL BIAS — HTS DOMINANCE")
    print(f"{'='*70}")

    total = len(df)
    meth_col = "Experiment_Methods"
    sys_col = "Experiment_System"

    # Records with any method annotation (not empty / "." / "n.a.")
    has_method = (
        df[meth_col].notna() &
        ~df[meth_col].astype(str).str.strip().isin(["", ".", "n.a."])
    )
    n_method_annotated = has_method.sum()

    # HTS records
    is_hts = df[meth_col].str.contains(
        r"(?i)high.?throughput|HTS|drug\s*screening", na=False
    )
    n_hts = is_hts.sum()

    # HTS within in vitro records that have methods
    is_in_vitro = df[sys_col].str.contains(r"(?i)in\s*vitro", na=False)
    in_vitro_with_methods = (is_in_vitro & has_method).sum()
    hts_in_vitro = (is_hts & is_in_vitro).sum()

    print(f"   Total records with method annotation:  {n_method_annotated:>6}  "
          f"({n_method_annotated/total*100:.1f}%)")
    print(f"   HTS records (total):                   {n_hts:>6}  "
          f"({n_hts/total*100:.1f}%)")
    print(f"   HTS as % of method-annotated:           {n_hts/n_method_annotated*100:>5.1f}%")
    print(f"")
    print(f"   In vitro records with methods:          {in_vitro_with_methods:>6}")
    print(f"   HTS within in vitro:                    {hts_in_vitro:>6}")
    print(f"   HTS as % of in-vitro-with-methods:      {hts_in_vitro/in_vitro_with_methods*100:>5.1f}%")

    # Detailed method value counts (top 10)
    print(f"\n   Top 10 Experiment_Methods:")
    meth_counts = df[meth_col].value_counts()
    for val, cnt in meth_counts.head(10).items():
        print(f"     {val}: {cnt} ({cnt/total*100:.1f}%)")

    return [
        {"metric": "Method-annotated records", "value": n_method_annotated,
         "pct": f"{n_method_annotated/total*100:.1f}%"},
        {"metric": "HTS records", "value": n_hts,
         "pct": f"{n_hts/total*100:.1f}%"},
        {"metric": "HTS % of method-annotated", "value": None,
         "pct": f"{n_hts/n_method_annotated*100:.1f}%"},
        {"metric": "HTS % of in-vitro-with-methods", "value": None,
         "pct": f"{hts_in_vitro/in_vitro_with_methods*100:.1f}%"},
    ]


# ---------------------------------------------------------------------------
# 4. Protein-level annotation coverage
# ---------------------------------------------------------------------------
def stat_protein_annotations(df):
    """Count records with protein-level annotations (Gene_Info column)."""
    print(f"\n{'='*70}")
    print(f"4. PROTEIN-LEVEL ANNOTATION COVERAGE")
    print(f"{'='*70}")

    total = len(df)

    # Gene_Info is the consolidated protein/gene annotation column
    n_gene_info = df["Gene_Info"].notna().sum()
    # Related_protein is the raw protein name
    n_related_protein = df["Related_protein"].notna().sum()
    # UniprotAC for UniProt accession
    n_uniprot = df["UniprotAC"].notna().sum()
    # FASTA sequence
    n_fasta = df["FASTA"].notna().sum() if "FASTA" in df.columns else 0
    # Function description
    n_function = df["Function"].notna().sum() if "Function" in df.columns else 0

    print(f"   Records with Gene_Info:             {n_gene_info:>6}  "
          f"({n_gene_info/total*100:.1f}%)")
    print(f"   Records with Related_protein:       {n_related_protein:>6}  "
          f"({n_related_protein/total*100:.1f}%)")
    print(f"   Records with UniprotAC:             {n_uniprot:>6}  "
          f"({n_uniprot/total*100:.1f}%)")
    print(f"   Records with FASTA:                 {n_fasta:>6}  "
          f"({n_fasta/total*100:.1f}%)" if total else "")
    print(f"   Records with Function:              {n_function:>6}  "
          f"({n_function/total*100:.1f}%)" if total else "")

    return {"metric": "Protein-level annotations (Gene_Info)",
            "value": n_gene_info, "pct": f"{n_gene_info/total*100:.1f}%"}


# ---------------------------------------------------------------------------
# 5. EC classification coverage
# ---------------------------------------------------------------------------
def stat_ec_classification(df):
    """Count records with Enzyme Commission classification."""
    print(f"\n{'='*70}")
    print(f"5. EC CLASSIFICATION COVERAGE")
    print(f"{'='*70}")

    total = len(df)
    n_ec = df["EC_Classification"].notna().sum()

    # Parse EC numbers for unique enzymes
    ec_values = df.loc[df["EC_Classification"].notna(), "EC_Classification"]
    unique_ec_numbers = set()
    for val in ec_values:
        # EC_Classification format: "2.7.1.95 | Transferase → Kinase → ..."
        ec_num = str(val).split(" | ")[0].strip()
        unique_ec_numbers.add(ec_num)

    print(f"   Records with EC_Classification:     {n_ec:>6}  "
          f"({n_ec/total*100:.1f}%)")
    print(f"   Unique EC numbers:                  {len(unique_ec_numbers):>6}")

    # Show EC number distribution
    print(f"\n   Top EC numbers by record count:")
    ec_record_counts = {}
    for val in ec_values:
        ec_num = str(val).split(" | ")[0].strip()
        ec_record_counts[ec_num] = ec_record_counts.get(ec_num, 0) + 1
    sorted_ecs = sorted(ec_record_counts.items(), key=lambda x: -x[1])
    for ec_num, cnt in sorted_ecs[:10]:
        print(f"     {ec_num}: {cnt}")

    return {"metric": "EC classification records", "value": n_ec,
            "pct": f"{n_ec/total*100:.1f}%"}


# ---------------------------------------------------------------------------
# 6–7. Michaelis constant (Km) coverage
# ---------------------------------------------------------------------------
def stat_km_coverage(df):
    """Count records with Michaelis constant and Km + combination data."""
    print(f"\n{'='*70}")
    print(f"6–7. MICHAELIS CONSTANT (Km) COVERAGE")
    print(f"{'='*70}")

    total = len(df)
    n_km = df["KM"].notna().sum()
    n_km_comb = df["KM_combination"].notna().sum()
    n_both = (df["KM"].notna() & df["KM_combination"].notna()).sum()

    print(f"   Records with KM:                    {n_km:>6}  "
          f"({n_km/total*100:.2f}%)")
    print(f"   Records with KM_combination:        {n_km_comb:>6}  "
          f"({n_km_comb/total*100:.2f}%)")
    print(f"   Records with both KM + KM_comb:     {n_both:>6}  "
          f"({n_both/total*100:.2f}%)")

    return [
        {"metric": "KM records", "value": n_km,
         "pct": f"{n_km/total*100:.2f}%"},
        {"metric": "KM_combination records", "value": n_km_comb,
         "pct": f"{n_km_comb/total*100:.2f}%"},
    ]


# ---------------------------------------------------------------------------
# 8. Unique annotated enzymes with kinetic data
# ---------------------------------------------------------------------------
def stat_unique_enzymes_with_km(df):
    """Count unique enzymes (by protein, gene, EC) that have Km data."""
    print(f"\n{'='*70}")
    print(f"8. UNIQUE ANNOTATED ENZYMES WITH KINETIC DATA")
    print(f"{'='*70}")

    km_df = df[df["KM"].notna()]

    n_unique_proteins = km_df["Related_protein"].nunique()
    n_unique_genes = km_df["Gene_Info"].nunique()
    n_unique_uniprot = km_df["UniprotAC"].nunique()

    # Unique EC numbers among Km records
    unique_ec_km = set()
    for val in km_df["EC_Classification"].dropna():
        ec_num = str(val).split(" | ")[0].strip()
        unique_ec_km.add(ec_num)

    print(f"   Unique proteins (Related_protein):  {n_unique_proteins}")
    print(f"   Unique gene annotations (Gene_Info): {n_unique_genes}")
    print(f"   Unique UniProt accessions:           {n_unique_uniprot}")
    print(f"   Unique EC numbers:                   {len(unique_ec_km)}")

    # List the enzymes
    print(f"\n   Enzymes with Km data (by Related_protein):")
    for prot in km_df["Related_protein"].dropna().unique():
        count = (km_df["Related_protein"] == prot).sum()
        print(f"     {prot}: {count} records")

    return {"metric": "Unique enzymes with Km (by protein)",
            "value": n_unique_proteins, "pct": None}


# ---------------------------------------------------------------------------
# 9. Unique drugs with kinetic data
# ---------------------------------------------------------------------------
def stat_drugs_with_km(df):
    """Count unique drugs total and those with kinetic data."""
    print(f"\n{'='*70}")
    print(f"9. UNIQUE DRUGS WITH KINETIC DATA")
    print(f"{'='*70}")

    total_drugs = df["Substance_Name"].nunique()
    km_drugs = df.loc[df["KM"].notna(), "Substance_Name"].nunique()

    print(f"   Total unique drugs:                 {total_drugs}")
    print(f"   Unique drugs with Km data:          {km_drugs}  "
          f"({km_drugs/total_drugs*100:.1f}%)")

    # Also count unique CIDs
    total_cids = df["PubChem_CID"].nunique()
    km_cids = df.loc[df["KM"].notna(), "PubChem_CID"].nunique()
    print(f"   Total unique PubChem CIDs:          {total_cids}")
    print(f"   Unique CIDs with Km data:           {km_cids}")

    return {"metric": "Drugs with Km / total drugs",
            "value": km_drugs, "pct": f"{km_drugs}/{total_drugs} ({km_drugs/total_drugs*100:.1f}%)"}


# ---------------------------------------------------------------------------
# 10. Source database distribution
# ---------------------------------------------------------------------------
def stat_source_distribution(df):
    """Distribution of records by source database."""
    print(f"\n{'='*70}")
    print(f"10. SOURCE DATABASE DISTRIBUTION")
    print(f"{'='*70}")

    total = len(df)
    if "Source_Dataset" not in df.columns:
        print("   Source_Dataset column not present in this file.")
        print("   (Available in intermediate pipeline files)")
        return [{"metric": "Source distribution", "value": None, "pct": "N/A"}]

    src_counts = df["Source_Dataset"].value_counts()
    for src, cnt in src_counts.items():
        print(f"   {src}: {cnt} ({cnt/total*100:.1f}%)")

    return [{"metric": f"Source: {src}", "value": cnt,
             "pct": f"{cnt/total*100:.1f}%"} for src, cnt in src_counts.items()]


# ---------------------------------------------------------------------------
# 11. Microbe diversity
# ---------------------------------------------------------------------------
def stat_microbe_diversity(df):
    """Microbe taxonomy diversity statistics."""
    print(f"\n{'='*70}")
    print(f"11. MICROBE DIVERSITY")
    print(f"{'='*70}")

    total = len(df)
    n_species = df["Microbe_Name"].nunique()
    n_genera = df["genus_name"].nunique() if "genus_name" in df.columns else 0
    n_families = df["family_name"].nunique() if "family_name" in df.columns else 0

    print(f"   Unique species:                     {n_species}")
    print(f"   Unique genera:                      {n_genera}")
    print(f"   Unique families:                    {n_families}")

    # Top 10 microbes by record count
    print(f"\n   Top 10 microbes by record count:")
    microbe_counts = df["Microbe_Name"].value_counts()
    for name, cnt in microbe_counts.head(10).items():
        print(f"     {name}: {cnt} ({cnt/total*100:.1f}%)")

    return {"metric": "Unique microbes", "value": n_species, "pct": None}


# ---------------------------------------------------------------------------
# 12. Drug diversity
# ---------------------------------------------------------------------------
def stat_drug_diversity(df):
    """Drug diversity statistics."""
    print(f"\n{'='*70}")
    print(f"12. DRUG DIVERSITY")
    print(f"{'='*70}")

    total = len(df)
    n_drugs = df["Substance_Name"].nunique()
    n_cids = df["PubChem_CID"].nunique()
    n_inchikeys = df["InChIKey"].nunique() if "InChIKey" in df.columns else 0

    print(f"   Unique drug names:                  {n_drugs}")
    print(f"   Unique PubChem CIDs:                {n_cids}")
    print(f"   Unique InChIKeys:                   {n_inchikeys}")

    # Substance category distribution
    if "Substance-Category" in df.columns:
        print(f"\n   Substance-Category distribution:")
        cat_counts = df["Substance-Category"].value_counts()
        for cat, cnt in cat_counts.items():
            print(f"     {cat}: {cnt} ({cnt/total*100:.1f}%)")

    # Top 10 drugs by record count
    print(f"\n   Top 10 drugs by record count:")
    drug_counts = df["Substance_Name"].value_counts()
    for name, cnt in drug_counts.head(10).items():
        print(f"     {name}: {cnt} ({cnt/total*100:.1f}%)")

    return {"metric": "Unique drugs", "value": n_drugs, "pct": None}


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def run_analysis():
    """Run all statistical analyses and produce a summary report."""
    print("=" * 70)
    print("STATISTICAL ANALYSIS OF MERGED MASI–MDIPID DATASET")
    print("=" * 70)

    df = load_dataset()

    results = []

    def add_result(r):
        """Add result(s) — handle both single dict and list of dicts."""
        if isinstance(r, list):
            results.extend(r)
        else:
            results.append(r)

    # Run all analyses
    add_result(stat_total_records(df))
    add_result(stat_experiment_system(df))
    add_result(stat_hts_bias(df))
    add_result(stat_protein_annotations(df))
    add_result(stat_ec_classification(df))
    add_result(stat_km_coverage(df))
    add_result(stat_unique_enzymes_with_km(df))
    add_result(stat_drugs_with_km(df))
    add_result(stat_source_distribution(df))
    add_result(stat_microbe_diversity(df))
    add_result(stat_drug_diversity(df))

    # Save summary CSV
    summary_df = pd.DataFrame(results)
    output_path = os.path.join(OUTPUT_DIR, "statistical_analysis_summary.csv")
    summary_df.to_csv(output_path, index=False)
    print(f"\n{'='*70}")
    print(f"SUMMARY TABLE SAVED: {output_path}")
    print(f"{'='*70}")

    # Print key findings matching manuscript
    print(f"\n{'='*70}")
    print(f"KEY FINDINGS (matching manuscript)")
    print(f"{'='*70}")
    total = len(df)
    n_in_vitro = df["Experiment_System"].str.contains(r"(?i)in\s*vitro", na=False).sum()
    is_in_vitro = df["Experiment_System"].str.contains(r"(?i)in\s*vitro", na=False)
    has_method = (
        df["Experiment_Methods"].notna() &
        ~df["Experiment_Methods"].astype(str).str.strip().isin(["", ".", "n.a."])
    )
    is_hts = df["Experiment_Methods"].str.contains(
        r"(?i)high.?throughput|HTS|drug\s*screening", na=False
    )
    in_vitro_with_methods = (is_in_vitro & has_method).sum()
    hts_in_vitro = (is_hts & is_in_vitro).sum()
    n_gene_info = df["Gene_Info"].notna().sum()
    n_ec = df["EC_Classification"].notna().sum()
    n_km = df["KM"].notna().sum()
    n_km_comb = df["KM_combination"].notna().sum()
    km_df = df[df["KM"].notna()]
    n_enzymes = km_df["Related_protein"].nunique()
    n_drugs_total = df["Substance_Name"].nunique()
    n_drugs_km = km_df["Substance_Name"].nunique()

    print(f"  Total records:                        {total:,}")
    print(f"  In vitro:                             {n_in_vitro} ({n_in_vitro/total*100:.1f}%)")
    print(f"  HTS % of in-vitro-with-methods:       {hts_in_vitro/in_vitro_with_methods*100:.1f}%")
    print(f"  Protein annotations (Gene_Info):      {n_gene_info} ({n_gene_info/total*100:.1f}%)")
    print(f"  EC classification:                    {n_ec} ({n_ec/total*100:.1f}%)")
    print(f"  Km values:                            {n_km} ({n_km/total*100:.2f}%)")
    print(f"  Km + combination:                     {n_km_comb} ({n_km_comb/total*100:.2f}%)")
    print(f"  Unique enzymes with Km:               {n_enzymes}")
    print(f"  Drugs with Km / total drugs:          {n_drugs_km}/{n_drugs_total} ({n_drugs_km/n_drugs_total*100:.1f}%)")

    return df


if __name__ == "__main__":
    run_analysis()
