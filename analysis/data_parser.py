"""
Data parser for DESeq2 results files.
Handles GFF-annotated CSV/TSV files with DESeq2 output columns.
"""

import pandas as pd
import numpy as np
import re
import io


def extract_attribute(attr_str, key):
    """Extract a value for a given key from a GFF-style Attributes string."""
    if pd.isna(attr_str):
        return None
    match = re.search(rf'{key}=([^;]+)', str(attr_str))
    return match.group(1) if match else None


def parse_deseq2_results(file_input, sep='\t'):
    """
    Parse DESeq2 results file and extract relevant columns.

    Parameters
    ----------
    file_input : str or file-like
        Path to the file or uploaded file buffer.
    sep : str
        Column separator (default: tab).

    Returns
    -------
    pd.DataFrame
        Parsed DataFrame with extracted gene information.
    """
    if isinstance(file_input, str):
        df = pd.read_csv(file_input, sep=sep)
    else:
        content = file_input.read()
        if isinstance(content, bytes):
            content = content.decode('utf-8')
        # Auto-detect separator
        first_line = content.split('\n')[0]
        if '\t' in first_line:
            sep = '\t'
        elif ',' in first_line and first_line.count(',') > 5:
            sep = ','
        df = pd.read_csv(io.StringIO(content), sep=sep)

    # Extract gene information from Attributes column
    if 'Attributes' in df.columns:
        df['locus_tag'] = df['Attributes'].apply(lambda x: extract_attribute(x, 'locus_tag'))
        df['gene_name'] = df['Attributes'].apply(lambda x: extract_attribute(x, 'gene'))
        df['product'] = df['Attributes'].apply(lambda x: extract_attribute(x, 'product'))
        df['protein_id'] = df['Attributes'].apply(lambda x: extract_attribute(x, 'protein_id'))

        # Extract GO terms
        def extract_go_terms(attr):
            if pd.isna(attr):
                return []
            match = re.search(r'Ontology_term=([^;]+)', str(attr))
            if match:
                return match.group(1).split(',')
            return []

        df['go_terms'] = df['Attributes'].apply(extract_go_terms)

    # Normalize column names for log2FoldChange and padj
    col_mapping = {}
    for col in df.columns:
        lower = col.lower().strip()
        if lower == 'log2foldchange' or lower == 'log2fc':
            col_mapping[col] = 'log2FoldChange'
        elif lower in ('padj', 'p.adj', 'pvalue_adjusted', 'fdr', 'qvalue'):
            col_mapping[col] = 'padj'
        elif lower == 'pvalue' or lower == 'p.value':
            col_mapping[col] = 'pvalue'
        elif lower == 'basemean':
            col_mapping[col] = 'baseMean'

    if col_mapping:
        df = df.rename(columns=col_mapping)

    # Convert numeric columns
    for col in ['log2FoldChange', 'padj', 'pvalue', 'baseMean']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Create a gene_id column (prefer locus_tag, fall back to gene_name, then index)
    if 'locus_tag' in df.columns:
        df['gene_id'] = df['locus_tag']
    elif 'gene_name' in df.columns:
        df['gene_id'] = df['gene_name']
    else:
        df['gene_id'] = df.index.astype(str)

    # Remove rows without gene_id
    df = df.dropna(subset=['gene_id'])
    # Remove duplicate gene IDs (keep first occurrence)
    df = df.drop_duplicates(subset=['gene_id'], keep='first')

    return df


def get_de_genes(df, padj_cutoff=0.05, log2fc_cutoff=1.0):
    """
    Get differentially expressed genes based on cutoffs.

    Returns
    -------
    tuple of (up_genes, down_genes, all_de_genes) as DataFrames
    """
    valid = df.dropna(subset=['log2FoldChange', 'padj'])

    up = valid[(valid['padj'] < padj_cutoff) & (valid['log2FoldChange'] > log2fc_cutoff)]
    down = valid[(valid['padj'] < padj_cutoff) & (valid['log2FoldChange'] < -log2fc_cutoff)]
    all_de = valid[(valid['padj'] < padj_cutoff) & (valid['log2FoldChange'].abs() > log2fc_cutoff)]

    return up, down, all_de


def get_gene_ranking(df, id_col='gene_id', fc_col='log2FoldChange'):
    """
    Create a gene ranking series for GSEA (sorted by log2FC, descending).

    Returns
    -------
    pd.Series
        Gene ranking indexed by gene_id.
    """
    valid = df.dropna(subset=[fc_col]).copy()
    valid = valid.drop_duplicates(subset=[id_col], keep='first')
    ranking = valid.set_index(id_col)[fc_col].sort_values(ascending=False)
    return ranking
