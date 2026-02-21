"""
COG (Clusters of Orthologous Groups) Analysis module.
Performs COG category enrichment analysis and functional classification.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import requests
import time
import re
import streamlit as st


# COG functional category definitions
COG_CATEGORIES = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}

# COG super-categories for grouping
COG_SUPERCATEGORIES = {
    "INFORMATION STORAGE AND PROCESSING": ["A", "B", "J", "K", "L"],
    "CELLULAR PROCESSES AND SIGNALING": ["D", "M", "N", "O", "T", "U", "V", "W", "Y", "Z"],
    "METABOLISM": ["C", "E", "F", "G", "H", "I", "P", "Q"],
    "POORLY CHARACTERIZED": ["R", "S", "X"],
}

# Color palette for COG categories
COG_COLORS = {
    "A": "#e6194b", "B": "#3cb44b", "C": "#ffe119", "D": "#4363d8",
    "E": "#f58231", "F": "#911eb4", "G": "#46f0f0", "H": "#f032e6",
    "I": "#bcf60c", "J": "#fabebe", "K": "#008080", "L": "#e6beff",
    "M": "#9a6324", "N": "#fffac8", "O": "#800000", "P": "#aaffc3",
    "Q": "#808000", "R": "#ffd8b1", "S": "#000075", "T": "#808080",
    "U": "#000000", "V": "#006400", "W": "#d7f9fa", "X": "#8c7088",
    "Y": "#c0c0c0", "Z": "#ff6347",
}


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_cog_from_kegg(org_code):
    """
    Attempt to fetch COG assignments from KEGG REST API.

    Returns
    -------
    dict
        {gene_id: [cog_ids]}
    """
    url = f"https://rest.kegg.jp/link/cog/{org_code}"
    try:
        time.sleep(0.35)
        resp = requests.get(url, timeout=60)
        if resp.status_code != 200:
            return {}
        links = {}
        for line in resp.text.strip().split('\n'):
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                gene_id = parts[0].replace(f'{org_code}:', '')
                cog_id = parts[1].replace('cog:', '')
                if gene_id not in links:
                    links[gene_id] = []
                links[gene_id].append(cog_id)
        return links
    except Exception:
        return {}


def parse_cog_annotation_file(file_input):
    """
    Parse a COG annotation file.

    Supported formats:
    - Tab-separated: gene_id <tab> COG_category (single letter or multiple)
    - eggNOG-mapper output format (columns with COG_category)
    - Simple two-column format

    Returns
    -------
    dict
        {gene_id: [cog_category_letters]}
    """
    import io

    if isinstance(file_input, str):
        df = pd.read_csv(file_input, sep='\t', comment='#')
    else:
        content = file_input.read()
        if isinstance(content, bytes):
            content = content.decode('utf-8')
        # Skip comment lines
        lines = [l for l in content.split('\n') if not l.startswith('#') and l.strip()]
        df = pd.read_csv(io.StringIO('\n'.join(lines)), sep='\t', header=None)

    mapping = {}

    if len(df.columns) >= 2:
        gene_col = df.columns[0]
        cog_col = df.columns[1]

        for _, row in df.iterrows():
            gene_id = str(row[gene_col]).strip()
            cog_str = str(row[cog_col]).strip()
            if cog_str and cog_str != 'nan' and cog_str != '-':
                # COG categories can be single letters or multiple (e.g., "COG0001" or "J" or "JK")
                cats = []
                for ch in cog_str:
                    if ch.upper() in COG_CATEGORIES:
                        cats.append(ch.upper())
                if cats:
                    mapping[gene_id] = cats

    return mapping


def infer_cog_from_products(df, product_col='product'):
    """
    Infer COG functional categories from gene product descriptions.
    This is a heuristic approach for when no COG annotation file is available.

    Returns
    -------
    dict
        {gene_id: [cog_category_letters]}
    """
    # Keyword-to-COG mapping (heuristic)
    keyword_mapping = {
        'C': ['dehydrogenase', 'oxidoreductase', 'cytochrome', 'electron transfer',
              'ATP synthase', 'NADH', 'succinate', 'fumarate', 'energy'],
        'E': ['amino acid', 'aminotransferase', 'transaminase', 'protease',
              'peptidase', 'amino acid transport', 'glutamate', 'aspartate',
              'threonine', 'lysine', 'methionine', 'serine', 'glycine', 'alanine',
              'glutamine', 'asparagine', 'proline', 'histidine', 'tryptophan',
              'tyrosine', 'phenylalanine', 'leucine', 'isoleucine', 'valine',
              'cysteine', 'arginine'],
        'F': ['nucleotide', 'purine', 'pyrimidine', 'kinase', 'nucleoside'],
        'G': ['carbohydrate', 'sugar', 'glycosyl', 'galactose', 'glucose',
              'mannose', 'fructose', 'xylose', 'PTS', 'phosphotransferase system'],
        'H': ['coenzyme', 'cofactor', 'biotin', 'thiamin', 'riboflavin',
              'folate', 'molybdopterin', 'cobalamin'],
        'I': ['lipid', 'fatty acid', 'acyl', 'lipase', 'phospholipid'],
        'J': ['ribosom', 'translation', 'tRNA', 'rRNA', 'aminoacyl',
              'elongation factor', 'initiation factor'],
        'K': ['transcription', 'transcriptional regulator', 'RNA polymerase',
              'sigma factor', 'repressor', 'activator', 'DNA-binding transcription'],
        'L': ['DNA repair', 'recombinase', 'helicase', 'ligase', 'replication',
              'DNA polymerase', 'topoisomerase', 'gyrase', 'recombination'],
        'M': ['cell wall', 'membrane', 'envelope', 'lipopolysaccharide', 'peptidoglycan',
              'murein', 'outer membrane', 'porin', 'LPS'],
        'N': ['flagell', 'motility', 'chemotaxis', 'pilus', 'fimbri'],
        'O': ['chaperone', 'protease', 'heat shock', 'protein folding',
              'ubiquitin', 'DnaK', 'DnaJ', 'GroEL', 'GroES', 'ClpB', 'proteasome'],
        'P': ['iron', 'zinc', 'manganese', 'copper', 'magnesium', 'phosphate',
              'sulfate', 'ion transport', 'siderophore', 'ABC transporter'],
        'T': ['signal transduction', 'sensor', 'two-component', 'histidine kinase',
              'response regulator', 'sensor kinase'],
        'U': ['secretion', 'type III', 'type IV', 'type II', 'type VI',
              'Sec', 'Tat', 'vesicular', 'export'],
        'V': ['defense', 'restriction', 'resistance', 'antimicrobial',
              'toxin-antitoxin', 'CRISPR'],
        'X': ['transposase', 'integrase', 'phage', 'prophage', 'insertion element',
              'mobile element', 'IS element'],
        'D': ['cell division', 'FtsZ', 'chromosome partition', 'septum',
              'cell cycle'],
        'Q': ['secondary metabolite', 'polyketide', 'nonribosomal peptide'],
    }

    mapping = {}
    if product_col not in df.columns:
        return mapping

    for _, row in df.iterrows():
        gene_id = str(row.get('gene_id', ''))
        product = str(row.get(product_col, '')).lower()
        if not product or product == 'nan':
            continue

        cats = set()
        for cog_cat, keywords in keyword_mapping.items():
            for kw in keywords:
                if kw.lower() in product:
                    cats.add(cog_cat)
                    break

        if not cats:
            # Default to S (function unknown) if product exists but no match
            cats.add('S')

        if gene_id and gene_id != 'nan':
            mapping[gene_id] = list(cats)

    return mapping


def run_cog_enrichment(de_gene_ids, all_gene_ids, cog_mapping, min_count=1):
    """
    Perform enrichment analysis on COG functional categories.

    Parameters
    ----------
    de_gene_ids : list
        Differentially expressed gene IDs.
    all_gene_ids : list
        All gene IDs (background).
    cog_mapping : dict
        {gene_id: [cog_category_letters]}

    Returns
    -------
    pd.DataFrame
        Enrichment results per COG category.
    """
    de_set = set(de_gene_ids)
    all_set = set(all_gene_ids)

    # Count genes per COG category in DE and background
    de_cog_counts = {}
    all_cog_counts = {}

    for gene in all_set:
        if gene in cog_mapping:
            for cat in cog_mapping[gene]:
                all_cog_counts[cat] = all_cog_counts.get(cat, 0) + 1

    for gene in de_set:
        if gene in cog_mapping:
            for cat in cog_mapping[gene]:
                de_cog_counts[cat] = de_cog_counts.get(cat, 0) + 1

    N = len(all_set)
    n = len(de_set)
    results = []

    for cat in sorted(COG_CATEGORIES.keys()):
        K = all_cog_counts.get(cat, 0)
        k = de_cog_counts.get(cat, 0)

        if K < min_count:
            continue

        # Fisher's exact test
        table = [[k, n - k], [K - k, max(0, N - K - n + k)]]
        try:
            odds_ratio, pvalue = stats.fisher_exact(table, alternative='greater')
        except ValueError:
            continue

        results.append({
            'COG_Category': cat,
            'Description': COG_CATEGORIES.get(cat, 'Unknown'),
            'DE_Count': k,
            'Background_Count': K,
            'DE_Total': n,
            'Background_Total': N,
            'Ratio_DE': round(k / n * 100, 2) if n > 0 else 0,
            'Ratio_BG': round(K / N * 100, 2) if N > 0 else 0,
            'Fold_Enrichment': round((k / n) / (K / N), 4) if K > 0 and n > 0 and N > 0 else 0,
            'pvalue': pvalue,
            'OddsRatio': round(odds_ratio, 4) if not np.isinf(odds_ratio) else 999.0,
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    _, padj, _, _ = multipletests(df['pvalue'].values, method='fdr_bh')
    df['p.adjust'] = padj
    df = df.sort_values('pvalue').reset_index(drop=True)
    return df


def get_cog_distribution(gene_ids, cog_mapping, label="All"):
    """
    Get the distribution of genes across COG categories.

    Returns
    -------
    pd.DataFrame
        Distribution table with COG category, count, percentage.
    """
    counts = {}
    mapped_count = 0

    for gene in gene_ids:
        if gene in cog_mapping:
            mapped_count += 1
            for cat in cog_mapping[gene]:
                counts[cat] = counts.get(cat, 0) + 1

    rows = []
    total = sum(counts.values())
    for cat in sorted(COG_CATEGORIES.keys()):
        count_val = counts.get(cat, 0)
        rows.append({
            'COG_Category': cat,
            'Description': COG_CATEGORIES[cat],
            'Count': count_val,
            'Percentage': round(count_val / total * 100, 2) if total > 0 else 0,
            'Group': label
        })

    return pd.DataFrame(rows)
