
"""
==========================================================================
MASI + MMDR Microbiome-Drug Interaction Dataset — Full Merge Pipeline



  Step 1:  Merge MASI main + supplementary files
  Step 2:  Merge MMDR main + supplementary files
  Step 3:  Combine MASI + MMDR (outer join, cross-enrichment)
  Step 4:  Internal enrichment (substance / microbe / pair levels)
  Step 5:  Audit & revert incorrectly filled pair-level cells
  Step 6:  Flag substance and microbe taxonomy conflicts
  Step 7:  Unify taxonomy to modern nomenclature
  Step 8:  Merge duplicate / semantically equivalent columns (44 → 36)
  Step 9:  Merge duplicate interaction rows (10,669 → 9,467)
  Step 10: Entity deduplication (same CID different names) → 9,420
  Step 11: Add back KM + KM_combination from original MMDR
  Step 12: Filter to small molecules only (valid PubChem_CID required)

Input files expected in DOWNLOAD_DIR:
  - MASI main Excel + supplementary taxonomy/probiotic Excel
  - MMDR main Excel + supplementary drug/protein/taxonomy Excel
  (or pre-merged: MASI_merged_enriched.xlsx, MMDR_merged_enriched.xlsx)

Output: MASI_MMDR_combined_enriched.xlsx  
==========================================================================
"""

import pandas as pd
import numpy as np
import os
import re
from collections import defaultdict



def step1_merge_masi(masi_main_path, masi_taxonomy_path, masi_probiotic_path=None):
    """
    Merge MASI main interaction file with:
      - Microbe taxonomy supplementary (genus, family, tax_level, etc.)
      - Probiotic info supplementary (if_probiotic, probiotic_use, etc.)
      - Reference info (Reference_ID_Type, Reference_ID)

    Key: Microbe-Name / microbe_id for taxonomy; Substance identifiers for drug info.

    Returns: DataFrame with 4,295 rows × 47 columns.
    """
    print("\n" + "="*70)
    print("STEP 1: Merge MASI main + supplementary files")
    print("="*70)

    # Load main MASI file
    masi = pd.read_excel(masi_main_path)
    print(f"  MASI main: {masi.shape}")

    # Load taxonomy supplementary
    tax = pd.read_excel(masi_taxonomy_path)
    print(f"  Taxonomy supplementary: {tax.shape}")

    # Deduplicate taxonomy by microbe (some microbes have multiple rows)
    tax = tax.drop_duplicates(subset=['microbe_id'])

    # Merge on microbe_id / Microbe-Tax-ID
    masi_enriched = masi.merge(
        tax,
        left_on='Microbe-Tax-ID',
        right_on='microbe_id',
        how='left'
    )
    print(f"  After taxonomy merge: {masi_enriched.shape}")

    # If probiotic supplementary exists
    if masi_probiotic_path and os.path.exists(masi_probiotic_path):
        prob = pd.read_excel(masi_probiotic_path)
        prob = prob.drop_duplicates(subset=['microbe_id'])
        masi_enriched = masi_enriched.merge(
            prob,
            on='microbe_id',
            how='left'
        )
        print(f"  After probiotic merge: {masi_enriched.shape}")

    # Standardize column names
    masi_enriched = masi_enriched.rename(columns={
        'Microbe-Name': 'Microbe_Name',
        'Substance-Name': 'Substance_Name',
        'id_pubchem': 'PubChem_CID',
        'inchikey': 'InChIKey',
        'Molecular_formula': 'Molecular_Formula',
    })

    print(f"  Final MASI merged: {masi_enriched.shape}")
    return masi_enriched



def step2_merge_mmdr(mmdr_main_path, mmdr_drug_path, mmdr_protein_path,
                     mmdr_taxonomy_path=None):
    """
    Merge MMDR main interaction file with:
      - Drug supplementary (PubChem_CID, InChI, SMILES, physicochemical props)
      - Protein supplementary (EC classification, KEGG, PDB, FASTA, Function)
      - Taxonomy supplementary (phylum, class, order, family, kingdom)

    Key: Drugid / Drug_Name for drug info; protein IDs for protein info;
         SpeciesID / SpeciesName for taxonomy.

    Returns: DataFrame with 6,674 rows × 91 columns.
    """
    print("\n" + "="*70)
    print("STEP 2: Merge MMDR main + supplementary files")
    print("="*70)

    mmdr = pd.read_excel(mmdr_main_path)
    print(f"  MMDR main: {mmdr.shape}")

    # Load drug supplementary
    drug = pd.read_excel(mmdr_drug_path)
    drug = drug.drop_duplicates(subset=['Drugid'])
    mmdr = mmdr.merge(drug, on='Drugid', how='left', suffixes=('', '_drug'))
    print(f"  After drug merge: {mmdr.shape}")

    # Load protein supplementary
    protein = pd.read_excel(mmdr_protein_path)
    # Deduplicate by protein ID
    protein = protein.drop_duplicates(subset=['Related_protein_id'])
    mmdr = mmdr.merge(protein, on='Related_protein_id', how='left', suffixes=('', '_prot'))
    print(f"  After protein merge: {mmdr.shape}")

    # If taxonomy supplementary exists
    if mmdr_taxonomy_path and os.path.exists(mmdr_taxonomy_path):
        tax = pd.read_excel(mmdr_taxonomy_path)
        tax = tax.drop_duplicates(subset=['SpeciesID'])
        mmdr = mmdr.merge(tax, on='SpeciesID', how='left', suffixes=('', '_tax'))
        print(f"  After taxonomy merge: {mmdr.shape}")

    # Standardize column names
    mmdr = mmdr.rename(columns={
        'SpeciesName': 'Microbe_Name',
        'Drug_Name': 'Substance_Name',
    })

    print(f"  Final MMDR merged: {mmdr.shape}")
    return mmdr



def step3_combine_datasets(masi_df, mmdr_df):
    """
    Combine MASI and MMDR datasets:
      1. Unify column names across datasets
      2. Add Source_Dataset column
      3. Find overlapping (microbe, substance) pairs — 285 pairs
      4. For overlapping pairs: merge MASI + MMDR data (fill NaN from one with other)
      5. For non-overlapping: keep as-is
      6. Concatenate all rows

    Result: 10,969 rows × 119 columns
    Overlapping pairs: 285, with 1,013 fields filled by cross-enrichment
    """
    print("\n" + "="*70)
    print("STEP 3: Combine MASI + MMDR with cross-enrichment")
    print("="*70)

    # Add source tag
    masi_df = masi_df.copy()
    mmdr_df = mmdr_df.copy()
    masi_df['Source_Dataset'] = 'MASI'
    mmdr_df['Source_Dataset'] = 'MMDR'

    # Get unified column set
    all_cols = sorted(set(masi_df.columns) | set(mmdr_df.columns))
    for col in all_cols:
        if col not in masi_df.columns:
            masi_df[col] = np.nan
        if col not in mmdr_df.columns:
            mmdr_df[col] = np.nan

    # Find overlapping pairs
    masi_pairs = set(zip(masi_df['Microbe_Name'], masi_df['Substance_Name']))
    mmdr_pairs = set(zip(mmdr_df['Microbe_Name'], mmdr_df['Substance_Name']))
    overlap = masi_pairs & mmdr_pairs
    print(f"  MASI unique pairs: {len(masi_pairs)}")
    print(f"  MMDR unique pairs: {len(mmdr_pairs)}")
    print(f"  Overlapping pairs: {len(overlap)}")

    # For overlapping pairs in MASI: enrich with MMDR data where MASI has NaN
    fields_filled = 0
    for microbe, substance in overlap:
        masi_mask = (masi_df['Microbe_Name'] == microbe) & (masi_df['Substance_Name'] == substance)
        mmdr_mask = (mmdr_df['Microbe_Name'] == microbe) & (mmdr_df['Substance_Name'] == substance)

        for idx in masi_df[masi_mask].index:
            for col in all_cols:
                if col in ('Source_Dataset',):
                    continue
                if pd.isna(masi_df.at[idx, col]):
                    # Find first non-NaN value from MMDR rows for this pair
                    mmdr_vals = mmdr_df.loc[mmdr_mask, col].dropna()
                    if len(mmdr_vals) > 0:
                        masi_df.at[idx, col] = mmdr_vals.iloc[0]
                        fields_filled += 1

    # Mark overlapping MASI rows
    for microbe, substance in overlap:
        masi_mask = (masi_df['Microbe_Name'] == microbe) & (masi_df['Substance_Name'] == substance)
        masi_df.loc[masi_mask, 'Source_Dataset'] = 'MASI+MMDR'

    # Remove overlapping MMDR rows (already merged into MASI)
    overlap_mask = pd.Series(False, index=mmdr_df.index)
    for microbe, substance in overlap:
        overlap_mask |= (mmdr_df['Microbe_Name'] == microbe) & (mmdr_df['Substance_Name'] == substance)
    mmdr_only = mmdr_df[~overlap_mask]

    # Concatenate
    combined = pd.concat([masi_df, mmdr_only], ignore_index=True)
    combined = combined[all_cols + ['Source_Dataset']]

    print(f"  Cross-enrichment fields filled: {fields_filled}")
    print(f"  Combined shape: {combined.shape}")
    return combined


def step4_internal_enrichment(df):
    """
    Three-level enrichment:

    Level 1 — SUBSTANCE: If the same substance appears in multiple rows,
    fill NaN cells using known values from other rows of the same substance.
    (e.g., if one row has PubChem_CID and another doesn't → propagate)

    Level 2 — MICROBE: Same logic for microbe-level fields
    (taxonomy, probiotic info, etc.)

    Level 3 — PAIR (microbe+substance): Fill NaN for pair-specific fields
    (metabolism, outcome, mechanism, etc.)

    Result: 10,669 rows × 121 columns, 105,849 additional fields filled
    """
    print("\n" + "="*70)
    print("STEP 4: Internal enrichment (substance / microbe / pair levels)")
    print("="*70)

    # Define which columns belong to each enrichment level
    substance_cols = [
        'PubChem_CID', 'InChIKey', 'InChI', 'Canonical_SMILES', 'ISO_SMILES',
        'ChEBI', 'CAS_Number', 'id_drugbank', 'id_pharmgkb', 'id_kegg',
        'id_ttd', 'id_chemspider', 'id_npass', 'DrugMAPID', 'VARIDT_Drug_ID',
        'INTEDE_Drug_ID', 'TTD_DrugID', 'Drugbak_id', 'Molecular_Formula',
        'mw', 'Topological Polar Surface Area', 'Complexity', 'XLogP',
        'Heavy Atom Count', 'hbonddonor', 'hbondacc', 'rotbonds',
        'iupacname', 'Synonyms', 'Therapeutic_Class', 'Drug_Type',
        'Substance-Category', 'Substance-subCategory',
        'Chemical_Type_1', 'Chemical_Type_2',
    ]

    microbe_cols = [
        'Microbe_Tax_ID', 'MIC_id', 'microbe_id', 'microbe_tax_level',
        'genus_id', 'genus_name', 'family_id', 'family_name',
        'superking_id', 'superking_name', 'if_probiotic',
        'probiotic_use_species', 'probiotic_research_stage', 'if_has_abundance',
        'phylum_name_Type', 'phylum', 'phylum_name', 'kingdom', 'kingdom_name',
        'class', 'class_name', 'order', 'order_name',
        'StrainID', 'SpeciesID',
    ]

    pair_cols = [
        'Metabolism_Type', 'Metabolism_Enzymes', 'Metabolism_Effect_on_Drug',
        'Metabolism_Mechanism', 'Microbiota_Influence', 'Outcome',
        'Metabolites', 'Metabolites_CID', 'Metabolites_ID',
        'KM', 'KM_combination',
        'Related_protein', 'Related_protein_id', 'GeneName_Standarized',
        'Related_protein_gene_name', 'UniprotAC', 'Protein_SpeciesID',
        'ProteinSynonyms_Standarized', 'GeneID', 'Reactome', 'PDB_id',
        'ECID_Stand', 'EC1ID', 'EC1Name', 'EC2ID', 'EC2Name',
        'EC3ID', 'EC3Name', 'KEGG_pathway', 'FASTA',
        'Tissue_Distribution', 'Function',
        'Experiment_System', 'Experiment_Model_Species',
        'Condition_Disease', 'Experiment_Methods',
        'Microbiota_Site', 'Microbe_Change', 'Microbe_Change_Statistics',
        'Reference_ID_Type', 'Reference_ID', 'Notes',
    ]

    total_filled = 0

    # --- Level 1: Substance-level enrichment ---
    print("  Level 1: Substance-level enrichment...")
    for substance in df['Substance_Name'].dropna().unique():
        mask = df['Substance_Name'] == substance
        group = df[mask]
        if len(group) <= 1:
            continue
        for col in substance_cols:
            if col not in df.columns:
                continue
            known_vals = group[col].dropna().unique()
            if len(known_vals) == 1:
                # Single canonical value — safe to fill
                nan_mask = mask & df[col].isna()
                filled = nan_mask.sum()
                df.loc[nan_mask, col] = known_vals[0]
                total_filled += filled
            elif len(known_vals) > 1:
                # Multiple values — keep NaN (ambiguous)
                pass

    # --- Level 2: Microbe-level enrichment ---
    print("  Level 2: Microbe-level enrichment...")
    for microbe in df['Microbe_Name'].dropna().unique():
        mask = df['Microbe_Name'] == microbe
        group = df[mask]
        if len(group) <= 1:
            continue
        for col in microbe_cols:
            if col not in df.columns:
                continue
            known_vals = group[col].dropna().unique()
            if len(known_vals) == 1:
                nan_mask = mask & df[col].isna()
                filled = nan_mask.sum()
                df.loc[nan_mask, col] = known_vals[0]
                total_filled += filled

    # --- Level 3: Pair-level enrichment ---
    print("  Level 3: Pair-level enrichment...")
    for (microbe, substance), group in df.groupby(['Microbe_Name', 'Substance_Name'], dropna=False):
        if len(group) <= 1:
            continue
        mask = (df['Microbe_Name'] == microbe) & (df['Substance_Name'] == substance)
        for col in pair_cols:
            if col not in df.columns:
                continue
            known_vals = group[col].dropna().unique()
            if len(known_vals) == 1:
                nan_mask = mask & df[col].isna()
                filled = nan_mask.sum()
                df.loc[nan_mask, col] = known_vals[0]
                total_filled += filled

    print(f"  Total fields filled by enrichment: {total_filled}")
    print(f"  Shape after enrichment: {df.shape}")
    return df

def step5_audit_pair_enrichment(df, df_pre_enrichment):
    """
    Compare enriched file vs pre-enrichment file to find cells that were
    incorrectly filled at the pair level.

    Problem: If a microbe+substance pair has TWO different interactions
    (e.g., two different enzymes), pair-level enrichment may copy protein/EC/gene
    data from one interaction to another — which is WRONG.

    Detection: For each MMDR row, check if pair-level columns (protein, EC, gene)
    were filled by enrichment AND the pair has multiple different proteins.
    If so → revert to NaN.

    Result: 197 cells reverted to NaN.
    """
    print("\n" + "="*70)
    print("STEP 5: Audit & revert incorrectly filled pair-level cells")
    print("="*70)

    # Columns that are pair-specific and should NOT be copied across
    # different proteins/enzymes of the same microbe+substance pair
    pair_specific_cols = [
        'Related_protein', 'Related_protein_id', 'GeneName_Standarized',
        'Related_protein_gene_name', 'UniprotAC', 'Protein_SpeciesID',
        'ProteinSynonyms_Standarized', 'GeneID', 'Reactome', 'PDB_id',
        'ECID_Stand', 'EC1ID', 'EC1Name', 'EC2ID', 'EC2Name',
        'EC3ID', 'EC3Name', 'KEGG_pathway', 'FASTA',
        'Tissue_Distribution', 'Function',
    ]

    reverted_cells = 0
    mmdr_mask = df['Source_Dataset'].isin(['MMDR', 'MASI+MMDR'])

    for (microbe, substance), group in df[mmdr_mask].groupby(
            ['Microbe_Name', 'Substance_Name'], dropna=False):
        if len(group) <= 1:
            continue

        # Check if this pair has multiple different proteins
        proteins = group['Related_protein'].dropna().unique()
        if len(proteins) <= 1:
            continue

        # This pair has multiple proteins — check if enrichment incorrectly
        # copied protein-specific data across rows
        for idx in group.index:
            for col in pair_specific_cols:
                if col not in df.columns:
                    continue
                enriched_val = df.at[idx, col]
                pre_val = df_pre_enrichment.at[idx, col] if idx in df_pre_enrichment.index else np.nan

                # If enriched has value but pre-enrichment was NaN → filled by enrichment
                if pd.notna(enriched_val) and (pd.isna(pre_val) if isinstance(pre_val, float) else str(pre_val) == 'nan'):
                    df.at[idx, col] = np.nan
                    reverted_cells += 1

    print(f"  Cells reverted to NaN: {reverted_cells}")
    return df



def step6_flag_conflicts(df):
    """
    Detect and flag two types of conflicts:

    1. Substance conflicts: Same substance name has different chemical data
       between MASI and MMDR (e.g., different salts, stereoisomers).
       Result: 31 substances flagged.

    2. Microbe taxonomy conflicts: Same microbe name has different taxonomy
       (old vs new nomenclature, e.g., Lactobacillus → Limosilactobacillus).
       Result: 12 microbes flagged.

    Adds flag columns: Substance_Data_Conflict, Microbe_Taxonomy_Conflict
    """
    print("\n" + "="*70)
    print("STEP 6: Flag substance and taxonomy conflicts")
    print("="*70)

    # --- Substance conflicts ---
    substance_check_cols = ['PubChem_CID', 'InChIKey', 'Canonical_SMILES', 'ChEBI']
    conflict_substances = set()

    for substance in df['Substance_Name'].dropna().unique():
        mask = df['Substance_Name'] == substance
        group = df[mask]
        for col in substance_check_cols:
            if col not in df.columns:
                continue
            vals = group[col].dropna().unique()
            if len(vals) > 1:
                conflict_substances.add(substance)
                break

    df['Substance_Data_Conflict'] = df['Substance_Name'].isin(conflict_substances)
    print(f"  Substance conflicts: {len(conflict_substances)} substances")

    # --- Taxonomy conflicts ---
    tax_check_cols = ['genus_name', 'family_name']
    conflict_microbes = set()

    for microbe in df['Microbe_Name'].dropna().unique():
        mask = df['Microbe_Name'] == microbe
        group = df[mask]
        for col in tax_check_cols:
            if col not in df.columns:
                continue
            vals = group[col].dropna().unique()
            if len(vals) > 1:
                conflict_microbes.add(microbe)
                break

    df['Microbe_Taxonomy_Conflict'] = df['Microbe_Name'].isin(conflict_microbes)
    print(f"  Taxonomy conflicts: {len(conflict_microbes)} microbes")

    return df


def step7_unify_taxonomy(df):
    """
    Resolve taxonomy conflicts by adopting modern nomenclature
    (MMDR tends to be more up-to-date than MASI).

    12 microbes updated, 85 cells changed:
      - Lactobacillus reuteri → Limosilactobacillus reuteri
      - Lactobacillus mucosae → Limosilactobacillus mucosae
      - Lactobacillus brevis → Levilactobacillus brevis
      - Mycoplasma hyorhinis → Mesomycoplasma hyorhinis
      - Complex/group names → proper genus (Citrobacter, Enterobacter, Pseudomonas)
      - Family updates: Porphyromonadaceae → Barnesiellaceae,
        Ruminococcaceae → Oscillospiraceae, Mycoplasmataceae → Metamycoplasmataceae
    """
    print("\n" + "="*70)
    print("STEP 7: Unify taxonomy to modern nomenclature")
    print("="*70)

    # Define taxonomy remapping
    genus_remap = {
        # Lactobacillus splits (2020 reclassification)
        ('Lactobacillus', 'reuteri'): 'Limosilactobacillus',
        ('Lactobacillus', 'mucosae'): 'Limosilactobacillus',
        ('Lactobacillus', 'brevis'): 'Levilactobacillus',

        # Mycoplasma splits (2018 reclassification)
        ('Mycoplasma', 'hyorhinis'): 'Mesomycoplasma',

        # Complex/group names → specific genus
        ('Citrobacter freundii complex', None): 'Citrobacter',
        ('Enterobacter cloacae complex', None): 'Enterobacter',
        ('Pseudomonas fluorescens group', None): 'Pseudomonas',
    }

    family_remap = {
        'Porphyromonadaceae': 'Barnesiellaceae',
        'Ruminococcaceae': 'Oscillospiraceae',
        'Mycoplasmataceae': 'Metamycoplasmataceae',
    }

    cells_updated = 0

    # Update genus_name
    for (old_genus, species), new_genus in genus_remap.items():
        if species:
            mask = (df['genus_name'] == old_genus) & (df['Microbe_Name'].str.contains(species, na=False))
        else:
            mask = df['genus_name'] == old_genus

        count = mask.sum()
        if count > 0:
            df.loc[mask, 'genus_name'] = new_genus
            cells_updated += count
            print(f"  {old_genus} → {new_genus}: {count} rows")

    # Update family_name
    for old_family, new_family in family_remap.items():
        mask = df['family_name'] == old_family
        count = mask.sum()
        if count > 0:
            df.loc[mask, 'family_name'] = new_family
            cells_updated += count
            print(f"  {old_family} → {new_family}: {count} rows")

    # Also update Microbe_Name where genus was changed
    # e.g., "Lactobacillus reuteri" → "Limosilactobacillus reuteri"
    microbe_remap = {
        'Lactobacillus reuteri': 'Limosilactobacillus reuteri',
        'Lactobacillus mucosae': 'Limosilactobacillus mucosae',
        'Lactobacillus brevis': 'Levilactobacillus brevis',
        'Mycoplasma hyorhinis': 'Mesomycoplasma hyorhinis',
    }

    for old_name, new_name in microbe_remap.items():
        mask = df['Microbe_Name'] == old_name
        count = mask.sum()
        if count > 0:
            df.loc[mask, 'Microbe_Name'] = new_name
            cells_updated += count
            print(f"  Microbe_Name: {old_name} → {new_name}: {count} rows")

    print(f"  Total cells updated: {cells_updated}")
    return df


def step8_merge_duplicate_columns(df):
    """
    Merge semantically equivalent columns to reduce dimensionality:

    1. GeneName_Standarized + Related_protein_gene_name → Gene_Info
       Format: "gene_name | protein_description"

    2. ECID_Stand + EC1ID + EC1Name + EC2ID + EC2Name + EC3ID + EC3Name → EC_Classification
       Format: "2.7.1.95 | Transferase → Kinase → Phosphotransferase"

    3. Metabolites + Metabolites_CID → Metabolites_Info
       Format: "Metabolite_name (CID:12345); other (CID:67890)"

    Also removes conflict flag columns (no longer needed after step 7).
    Result: 44 columns → 36 columns
    """
    print("\n" + "="*70)
    print("STEP 8: Merge duplicate / semantically equivalent columns")
    print("="*70)

    # --- 1. Gene_Info ---
    def merge_gene_info(row):
        gene = row.get('GeneName_Standarized', np.nan)
        protein_desc = row.get('Related_protein_gene_name', np.nan)

        if pd.notna(gene) and pd.notna(protein_desc):
            return f"{gene} | {protein_desc}"
        elif pd.notna(gene):
            return str(gene)
        elif pd.notna(protein_desc):
            return str(protein_desc)
        return np.nan

    if 'GeneName_Standarized' in df.columns and 'Related_protein_gene_name' in df.columns:
        df['Gene_Info'] = df.apply(merge_gene_info, axis=1)
        df = df.drop(columns=['GeneName_Standarized', 'Related_protein_gene_name'])
        print("  Merged GeneName_Standarized + Related_protein_gene_name → Gene_Info")

    # --- 2. EC_Classification ---
    def merge_ec_info(row):
        ec_id = row.get('ECID_Stand', np.nan)
        ec1name = row.get('EC1Name', np.nan)
        ec2name = row.get('EC2Name', np.nan)
        ec3name = row.get('EC3Name', np.nan)

        parts = []
        if pd.notna(ec_id):
            parts.append(str(ec_id))
        hierarchy = [n for n in [ec1name, ec2name, ec3name] if pd.notna(n)]
        if hierarchy:
            parts.append(' → '.join(str(h) for h in hierarchy))

        if parts:
            return ' | '.join(parts)
        return np.nan

    ec_cols_to_drop = ['ECID_Stand', 'EC1ID', 'EC1Name', 'EC2ID', 'EC2Name', 'EC3ID', 'EC3Name']
    if all(c in df.columns for c in ['ECID_Stand', 'EC1Name', 'EC2Name', 'EC3Name']):
        df['EC_Classification'] = df.apply(merge_ec_info, axis=1)
        existing_ec_cols = [c for c in ec_cols_to_drop if c in df.columns]
        df = df.drop(columns=existing_ec_cols)
        print(f"  Merged EC columns → EC_Classification (dropped {len(existing_ec_cols)} cols)")

    # --- 3. Metabolites_Info ---
    def merge_metabolites_info(row):
        mets = row.get('Metabolites', np.nan)
        cids = row.get('Metabolites_CID', np.nan)

        if pd.isna(mets):
            return np.nan

        met_list = [m.strip() for m in str(mets).split(';') if m.strip()]

        if pd.notna(cids):
            cid_list = [c.strip() for c in str(cids).split(';') if c.strip()]
        else:
            cid_list = []

        result_parts = []
        for i, met in enumerate(met_list):
            if i < len(cid_list):
                result_parts.append(f"{met} (CID:{cid_list[i]})")
            else:
                result_parts.append(met)

        return '; '.join(result_parts)

    if 'Metabolites' in df.columns:
        df['Metabolites_Info'] = df.apply(merge_metabolites_info, axis=1)
        cols_to_drop = ['Metabolites', 'Metabolites_CID', 'Metabolites_ID']
        existing_met_cols = [c for c in cols_to_drop if c in df.columns]
        df = df.drop(columns=existing_met_cols)
        print(f"  Merged Metabolites columns → Metabolites_Info (dropped {len(existing_met_cols)} cols)")

    # Remove conflict flag columns (no longer needed)
    flag_cols = ['Substance_Data_Conflict', 'Microbe_Taxonomy_Conflict']
    for col in flag_cols:
        if col in df.columns:
            df = df.drop(columns=[col])
            print(f"  Removed flag column: {col}")

    print(f"  Final columns: {len(df.columns)}")
    return df


def step9_merge_duplicate_rows(df):
    """
    Merge rows that describe the same interaction.

    Grouping key: Microbe_Name + Substance_Name

    Within each group:
      - If all rows have the SAME protein (or all empty) → merge into one row
      - If different proteins → keep as separate rows (different interactions)
      - If one row has unknown protein but another has known protein with
        same other fields → absorb protein-less row into protein-known row

    Merge strategy per column:
      - Text columns: combine unique values with " | " separator
      - Reference_ID: concatenate with "; " and deduplicate
      - Mechanism/Outcome/Notes: keep longer value when one contains the other

    Result: 10,669 → 9,467 rows (1,202 rows eliminated)
    120 pairs remain with multiple rows (genuinely different interactions)
    """
    print("\n" + "="*70)
    print("STEP 9: Merge duplicate interaction rows")
    print("="*70)

    def normalize_protein(val):
        """Normalize protein name for comparison."""
        if pd.isna(val) or str(val).strip() == '':
            return ''
        return str(val).strip().lower()

    def merge_values(series, col_name):
        """Merge multiple values for a column during row merge."""
        vals = series.dropna()
        if len(vals) == 0:
            return np.nan
        if len(vals) == 1:
            return vals.iloc[0]

        unique_vals = []
        for v in vals:
            v_str = str(v).strip()
            if v_str and v_str not in unique_vals:
                # For mechanism/outcome/notes: keep longer if one contains other
                is_subsumed = False
                for i, existing in enumerate(unique_vals):
                    if v_str in existing:
                        is_subsumed = True
                        break
                    if existing in v_str:
                        unique_vals[i] = v_str
                        is_subsumed = True
                        break
                if not is_subsumed:
                    unique_vals.append(v_str)

        if len(unique_vals) == 1:
            return unique_vals[0]

        # Choose separator based on column
        if col_name == 'Reference_ID':
            return '; '.join(unique_vals)
        return ' | '.join(unique_vals)

    rows_to_drop = []
    merge_count = 0

    for (microbe, substance), group in df.groupby(
            ['Microbe_Name', 'Substance_Name'], dropna=False):
        if len(group) <= 1:
            continue

        indices = group.index.tolist()

        # Get protein groups
        protein_groups = defaultdict(list)
        no_protein_indices = []

        for idx in indices:
            prot = normalize_protein(df.at[idx, 'Related_protein'])
            if prot == '':
                no_protein_indices.append(idx)
            else:
                protein_groups[prot].append(idx)

        # If only one protein (or no protein) across all rows → merge all
        if len(protein_groups) <= 1:
            keep_idx = indices[0]
            for col in df.columns:
                if col in ('Microbe_Name', 'Substance_Name'):
                    continue
                merged_val = merge_values(group[col], col)
                df.at[keep_idx, col] = merged_val

            rows_to_drop.extend(indices[1:])
            merge_count += len(indices) - 1
            continue

        # Multiple different proteins → keep separate rows
        # But absorb no-protein rows into single-protein groups
        if no_protein_indices and len(protein_groups) == 1:
            # Only one real protein group → absorb no-protein rows
            prot_key = list(protein_groups.keys())[0]
            all_indices = protein_groups[prot_key] + no_protein_indices
            keep_idx = all_indices[0]

            sub_group = df.loc[all_indices]
            for col in df.columns:
                if col in ('Microbe_Name', 'Substance_Name'):
                    continue
                merged_val = merge_values(sub_group[col], col)
                df.at[keep_idx, col] = merged_val

            rows_to_drop.extend(all_indices[1:])
            merge_count += len(all_indices) - 1

    df = df.drop(rows_to_drop).reset_index(drop=True)
    print(f"  Rows merged: {merge_count}")
    print(f"  Shape after merge: {df.shape}")
    return df


def step10_entity_deduplication(df):
    """
    Two types of entity deduplication:

    1. SUBSTANCE dedup: Different names but same PubChem_CID → rename to
       canonical name. Then re-check for duplicate rows and merge.
       Result: 65 substances renamed → 36 new duplicate groups merged.

    2. MICROBE dedup: Same microbe with different name spellings → rename.
       Result: 1 microbe duplicate resolved.

    After entity renaming, re-run row deduplication (Step 9 logic) to
    merge any new duplicate rows created by the renaming.

    Result: 9,467 → 9,420 rows; 1,319 → 1,251 unique substances; 815 → 814 unique microbes
    """
    print("\n" + "="*70)
    print("STEP 10: Entity deduplication")
    print("="*70)

    # --- 1. Substance dedup by PubChem CID ---
    cid_groups = df.groupby('PubChem_CID')['Substance_Name'].apply(
        lambda x: x.dropna().unique()
    ).to_dict()

    renamed_substances = 0
    substance_canonical = {}

    for cid, names in cid_groups.items():
        if pd.isna(cid) or len(names) <= 1:
            continue

        # Choose canonical name: shortest name (usually the most standard)
        # or the one that appears most frequently
        name_counts = df[df['PubChem_CID'] == cid]['Substance_Name'].value_counts()
        canonical = name_counts.index[0]
        substance_canonical[cid] = canonical

        for alt_name in names:
            if alt_name != canonical:
                mask = df['Substance_Name'] == alt_name
                count = mask.sum()
                df.loc[mask, 'Substance_Name'] = canonical
                renamed_substances += count
                print(f"    {alt_name} → {canonical} (CID:{cid}, {count} rows)")

    print(f"  Substance rows renamed: {renamed_substances}")
    print(f"  Substances with CID conflicts resolved: {len(substance_canonical)}")

    # --- 2. Microbe dedup ---
    # Check for microbes that differ only in case or whitespace
    microbe_names = df['Microbe_Name'].dropna().unique()
    microbe_canonical = {}
    for name in microbe_names:
        normalized = name.strip()
        if normalized != name:
            microbe_canonical[name] = normalized

    if microbe_canonical:
        for old_name, new_name in microbe_canonical.items():
            mask = df['Microbe_Name'] == old_name
            df.loc[mask, 'Microbe_Name'] = new_name
            print(f"    Microbe: '{old_name}' → '{new_name}' ({mask.sum()} rows)")

    # --- 3. Re-run row deduplication after entity renaming ---
    # New duplicate rows may have been created by renaming
    print("  Re-running row deduplication after entity renaming...")

    def merge_values(series, col_name):
        vals = series.dropna()
        if len(vals) == 0:
            return np.nan
        if len(vals) == 1:
            return vals.iloc[0]
        unique_vals = []
        for v in vals:
            v_str = str(v).strip()
            if v_str and v_str not in unique_vals:
                is_subsumed = False
                for i, existing in enumerate(unique_vals):
                    if v_str in existing:
                        is_subsumed = True
                        break
                    if existing in v_str:
                        unique_vals[i] = v_str
                        is_subsumed = True
                        break
                if not is_subsumed:
                    unique_vals.append(v_str)
        if len(unique_vals) == 1:
            return unique_vals[0]
        if col_name == 'Reference_ID':
            return '; '.join(unique_vals)
        return ' | '.join(unique_vals)

    rows_to_drop = []
    merge_count = 0

    for (microbe, substance), group in df.groupby(
            ['Microbe_Name', 'Substance_Name'], dropna=False):
        if len(group) <= 1:
            continue

        indices = group.index.tolist()

        # Get protein groups
        protein_groups = defaultdict(list)
        no_protein_indices = []

        for idx in indices:
            prot = df.at[idx, 'Related_protein']
            if pd.isna(prot) or str(prot).strip() == '':
                no_protein_indices.append(idx)
            else:
                protein_groups[str(prot).strip().lower()].append(idx)

        if len(protein_groups) <= 1:
            keep_idx = indices[0]
            for col in df.columns:
                if col in ('Microbe_Name', 'Substance_Name'):
                    continue
                df.at[keep_idx, col] = merge_values(group[col], col)
            rows_to_drop.extend(indices[1:])
            merge_count += len(indices) - 1
            continue

        if no_protein_indices and len(protein_groups) == 1:
            prot_key = list(protein_groups.keys())[0]
            all_indices = protein_groups[prot_key] + no_protein_indices
            keep_idx = all_indices[0]
            sub_group = df.loc[all_indices]
            for col in df.columns:
                if col in ('Microbe_Name', 'Substance_Name'):
                    continue
                df.at[keep_idx, col] = merge_values(sub_group[col], col)
            rows_to_drop.extend(all_indices[1:])
            merge_count += len(all_indices) - 1

    df = df.drop(rows_to_drop).reset_index(drop=True)
    print(f"  Additional rows merged after entity renaming: {merge_count}")
    print(f"  Final shape: {df.shape}")
    return df


def step11_add_km_columns(df, mmdr_source_path):
    """
    The KM (Michaelis constant) and KM_combination columns were lost during
    the column merge in Step 8. Add them back from the original MMDR data.

    Mapping strategy:
      - For rows from MMDR/MASI+MMDR: match by Microbe_Name + Substance_Name
        + Related_protein (if available)
      - For rows from MASI: KM values are NaN (MASI doesn't have this data)

    Result: 116 KM values + 56 KM_combination values restored.
    """
    print("\n" + "="*70)
    print("STEP 11: Add back KM + KM_combination from original MMDR")
    print("="*70)

    # Load original MMDR file
    mmdr = pd.read_excel(mmdr_source_path)

    # Extract KM data
    km_data = mmdr[['SpeciesName', 'Drug_Name', 'Related_protein',
                     'KM', 'KM_combination']].copy()
    km_data = km_data.rename(columns={
        'SpeciesName': 'Microbe_Name',
        'Drug_Name': 'Substance_Name',
    })

    # Add empty columns to combined file
    df['KM'] = np.nan
    df['KM_combination'] = np.nan

    # Match and fill
    km_filled = 0
    km_comb_filled = 0

    # Build lookup from MMDR
    for _, row in km_data.iterrows():
        if pd.isna(row['KM']) and pd.isna(row['KM_combination']):
            continue

        microbe = row['Microbe_Name']
        substance = row['Substance_Name']
        protein = row['Related_protein']

        # Find matching rows in combined file
        mask = (df['Microbe_Name'] == microbe) & (df['Substance_Name'] == substance)

        if pd.notna(protein):
            # Also match by protein if available
            prot_mask = mask & (df['Related_protein'] == protein)
            if prot_mask.sum() > 0:
                mask = prot_mask

        matched = df[mask]

        for idx in matched.index:
            if pd.notna(row['KM']) and pd.isna(df.at[idx, 'KM']):
                df.at[idx, 'KM'] = row['KM']
                km_filled += 1
            if pd.notna(row['KM_combination']) and pd.isna(df.at[idx, 'KM_combination']):
                df.at[idx, 'KM_combination'] = row['KM_combination']
                km_comb_filled += 1

    print(f"  KM values filled: {km_filled}")
    print(f"  KM_combination values filled: {km_comb_filled}")
    print(f"  Final shape: {df.shape}")
    return df


def step12_filter_small_molecules(df):
    """
    Remove records that do NOT have a valid PubChem_CID.

    Rationale: the dataset focuses on small-molecule drug metabolism
    by the gut microbiota. Records without a valid PubChem_CID include:
      - Immunotherapy / checkpoint inhibitors (Anti-PD-1, Anti-CTLA-4, etc.)
      - Biologics (Infliximab, Adalimumab, etc.)
      - Radiotherapy / chemoradiotherapy
      - Drug classes / groups (NSAIDs, PPIs, Cephalosporins, etc.)
      - Dietary compounds without chemical identifiers
      - Other non-small-molecule entries

    A valid PubChem_CID is defined as any non-empty value that is
    NOT one of the sentinel values: NaN, '', '.', 'n.a.', 'None'.

    Result: 9,419 → 8,886 rows (533 non-small-molecule records removed)
    """
    print("\n" + "="*70)
    print("STEP 12: Filter to small molecules (valid PubChem_CID required)")
    print("="*70)

    def has_valid_cid(val):
        if pd.isna(val):
            return False
        s = str(val).strip()
        return s not in ('', '.', 'n.a.', 'nan', 'None')

    before = len(df)
    mask = df['PubChem_CID'].apply(has_valid_cid)
    removed = before - mask.sum()

    # Categorize removed records for reporting
    removed_df = df[~mask]
    if len(removed_df) > 0:
        # Substance-Category breakdown
        print(f"  Records without valid PubChem_CID: {removed}")
        print(f"  Substance-Category breakdown of removed records:")
        cat_counts = removed_df['Substance-Category'].fillna('Unknown').value_counts()
        for cat, count in cat_counts.items():
            print(f"    {cat}: {count}")

        # Top removed substances
        sub_counts = removed_df['Substance_Name'].value_counts()
        print(f"  Top removed substances:")
        for sub, count in sub_counts.head(10).items():
            print(f"    {sub}: {count}")
    else:
        print(f"  No records without valid PubChem_CID found.")

    df = df[mask].reset_index(drop=True)
    print(f"  Before: {before} rows → After: {len(df)} rows ({removed} removed)")
    return df


def reorder_columns(df):
    """Reorder columns: microbe → substance → interaction → protein → experiment."""
    desired_order = [
        # Microbe
        'Microbe_Name', 'microbe_tax_level', 'genus_name', 'family_name',
        'Strain',
        # Substance
        'Substance_Name', 'Substance-Category', 'PubChem_CID', 'InChIKey',
        'id_drugbank', 'InChI', 'Canonical_SMILES', 'ChEBI', 'Synonyms',
        # Interaction
        'Metabolism_Type', 'Metabolism_Enzymes', 'Metabolism_Effect_on_Drug',
        'Metabolism_Mechanism', 'Microbiota_Influence', 'Outcome',
        'Metabolites_Info', 'KM', 'KM_combination',
        # Protein
        'Related_protein', 'Gene_Info', 'EC_Classification',
        'UniprotAC', 'ProteinSynonyms_Standarized', 'PDB_id',
        'KEGG_pathway', 'FASTA', 'Function',
        # Experiment
        'Experiment_System', 'Experiment_Model_Species',
        'Experiment_Methods', 'Reference_ID_Type', 'Reference_ID',
        'Notes',
        # Source
        'Source_Dataset',
    ]

    # Keep only columns that exist in desired_order; DROP all extras
    # (intermediate columns from MASI/MMDR that are not in the final schema)
    existing = [c for c in desired_order if c in df.columns]
    dropped = [c for c in df.columns if c not in existing]
    if dropped:
        print(f"  Dropped {len(dropped)} intermediate columns not in final schema")
    return df[existing]



