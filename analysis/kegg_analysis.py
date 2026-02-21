"""
KEGG Pathway Enrichment Analysis module.
Provides ORA (Over-Representation Analysis) and GSEA (Gene Set Enrichment Analysis)
for KEGG pathways using the KEGG REST API.
"""

import requests
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import time
import streamlit as st


KEGG_BASE = "https://rest.kegg.jp"
REQUEST_DELAY = 0.35  # seconds between API calls to respect rate limits


def _kegg_get(endpoint, max_retries=3):
    """Make a KEGG REST API request with retries and rate limiting."""
    url = f"{KEGG_BASE}/{endpoint}"
    for attempt in range(max_retries):
        try:
            time.sleep(REQUEST_DELAY)
            resp = requests.get(url, timeout=60)
            if resp.status_code == 200:
                return resp.text
            elif resp.status_code == 403:
                time.sleep(5)
                continue
        except requests.exceptions.RequestException:
            time.sleep(2)
    return ""


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_kegg_gene_list(org_code):
    """
    Fetch all genes for an organism from KEGG.

    Returns
    -------
    dict
        {gene_id: description}
        Description is the last field of the KEGG list line, e.g.,
        "thrA; bifunctional aspartokinase I..." or "hypothetical protein"
    """
    text = _kegg_get(f"list/{org_code}")
    genes = {}
    for line in text.strip().split('\n'):
        if not line.strip():
            continue
        # KEGG list format: org:gene_id\tCDS\tstart..end\tdescription
        # or sometimes:     org:gene_id\tdescription
        all_parts = line.split('\t')
        if len(all_parts) < 2:
            continue
        gene_id = all_parts[0].replace(f'{org_code}:', '')
        # The description with gene symbol is always the LAST tab field
        description = all_parts[-1].strip()
        genes[gene_id] = description
    return genes


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_kegg_pathways(org_code):
    """
    Fetch all pathway definitions for an organism.

    Returns
    -------
    dict
        {pathway_id: pathway_name}
    """
    text = _kegg_get(f"list/pathway/{org_code}")
    pathways = {}
    for line in text.strip().split('\n'):
        if not line.strip():
            continue
        parts = line.split('\t', 1)
        if len(parts) >= 2:
            pw_id = parts[0].replace('path:', '')
            # Remove organism suffix from pathway name
            pw_name = parts[1].split(' - ')[0].strip()
            pathways[pw_id] = pw_name
    return pathways


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_kegg_gene_pathway_links(org_code):
    """
    Fetch gene-pathway association links.

    Returns
    -------
    dict
        {gene_id: [pathway_ids]}
    """
    text = _kegg_get(f"link/pathway/{org_code}")
    links = {}
    for line in text.strip().split('\n'):
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            gene_id = parts[0].replace(f'{org_code}:', '')
            pw_id = parts[1].replace('path:', '')
            if gene_id not in links:
                links[gene_id] = []
            links[gene_id].append(pw_id)
    return links


def build_kegg_id_mapping(input_gene_ids, kegg_genes, org_code,
                         gene_name_map=None, product_map=None):
    """
    Map input gene IDs to KEGG gene IDs using multiple strategies.

    Strategies (in order):
    1. Exact match on gene ID
    2. Case-insensitive match on gene ID
    3. Gene name/symbol match (input gene_name → KEGG gene symbol from description)
    4. Product description match (partial)
    5. Partial ID match (fallback)

    Parameters
    ----------
    input_gene_ids : list
        Gene IDs from the input data (e.g., locus_tags).
    kegg_genes : dict
        {kegg_gene_id: description} from KEGG API.
    org_code : str
        KEGG organism code.
    gene_name_map : dict, optional
        {input_gene_id: gene_name/symbol} for gene name-based matching.
    product_map : dict, optional
        {input_gene_id: product_description} for product-based matching.

    Returns
    -------
    dict
        {input_gene_id: kegg_gene_id}
    """
    mapping = {}
    kegg_ids = list(kegg_genes.keys())
    kegg_ids_lower = {k.lower(): k for k in kegg_ids}

    # Build KEGG gene symbol → KEGG gene ID lookup
    kegg_symbol_to_id = {}  # gene_symbol_lower → kegg_gene_id
    kegg_product_to_id = {}  # product_keyword_lower → kegg_gene_id
    for kid, desc in kegg_genes.items():
        if ';' in desc:
            gene_sym = desc.split(';')[0].strip()
            kegg_symbol_to_id[gene_sym.lower()] = kid
            # Also store product part
            product_part = desc.split(';', 1)[1].strip().lower()
            if len(product_part) > 5:  # meaningful product description
                kegg_product_to_id[product_part] = kid
        else:
            # No semicolon — entire desc is product
            desc_lower = desc.strip().lower()
            if len(desc_lower) > 5:
                kegg_product_to_id[desc_lower] = kid

    if gene_name_map is None:
        gene_name_map = {}
    if product_map is None:
        product_map = {}

    for gid in input_gene_ids:
        if gid is None or pd.isna(gid):
            continue
        gid_str = str(gid).strip()

        # Strategy 1: Exact match on gene ID
        if gid_str in kegg_genes:
            mapping[gid_str] = gid_str
            continue

        # Strategy 2: Case-insensitive match on gene ID
        if gid_str.lower() in kegg_ids_lower:
            mapping[gid_str] = kegg_ids_lower[gid_str.lower()]
            continue

        # Strategy 3: Gene name/symbol match (PRIMARY for cross-strain)
        gene_name = gene_name_map.get(gid_str)
        if gene_name and str(gene_name).lower() in kegg_symbol_to_id:
            mapping[gid_str] = kegg_symbol_to_id[str(gene_name).lower()]
            continue

        # Strategy 4: Input gene_id itself as a gene symbol
        if gid_str.lower() in kegg_symbol_to_id:
            mapping[gid_str] = kegg_symbol_to_id[gid_str.lower()]
            continue

        # Strategy 5: Product description match
        product = product_map.get(gid_str)
        if product:
            product_lower = str(product).strip().lower()
            if len(product_lower) > 10:  # meaningful product
                for kegg_prod, kid in kegg_product_to_id.items():
                    if product_lower == kegg_prod or (
                        len(product_lower) > 15 and product_lower in kegg_prod
                    ):
                        mapping[gid_str] = kid
                        break
                if gid_str in mapping:
                    continue

    return mapping


def run_kegg_ora(de_gene_ids, all_gene_ids, gene_pathway_links, pathways,
                 id_mapping=None, min_size=3, max_size=500):
    """
    Perform Over-Representation Analysis (ORA) for KEGG pathways.

    Parameters
    ----------
    de_gene_ids : list
        Differentially expressed gene IDs.
    all_gene_ids : list
        All gene IDs (background).
    gene_pathway_links : dict
        {gene_id: [pathway_ids]}
    pathways : dict
        {pathway_id: pathway_name}
    id_mapping : dict, optional
        {input_gene_id: kegg_gene_id}
    min_size : int
        Minimum pathway size.
    max_size : int
        Maximum pathway size.

    Returns
    -------
    pd.DataFrame
        ORA results with columns: Pathway_ID, Description, GeneRatio, BgRatio,
        pvalue, p.adjust, OddsRatio, Count, Genes
    """
    # Apply ID mapping if provided
    if id_mapping:
        de_mapped = {id_mapping[g] for g in de_gene_ids if g in id_mapping}
        all_mapped = {id_mapping[g] for g in all_gene_ids if g in id_mapping}
    else:
        de_mapped = set(de_gene_ids)
        all_mapped = set(all_gene_ids)

    # Build pathway -> gene set
    pathway_gene_sets = {}
    for gene, pw_list in gene_pathway_links.items():
        if gene in all_mapped:
            for pw in pw_list:
                if pw not in pathway_gene_sets:
                    pathway_gene_sets[pw] = set()
                pathway_gene_sets[pw].add(gene)

    N = len(all_mapped)
    n = len(de_mapped)
    results = []

    for pw_id, pw_genes in pathway_gene_sets.items():
        K = len(pw_genes)
        if K < min_size or K > max_size:
            continue
        k = len(pw_genes.intersection(de_mapped))
        if k == 0:
            continue

        # Fisher's exact test (one-sided, greater)
        table = [[k, n - k], [K - k, max(0, N - K - n + k)]]
        try:
            odds_ratio, pvalue = stats.fisher_exact(table, alternative='greater')
        except ValueError:
            continue

        pw_name = pathways.get(pw_id, pw_id)
        results.append({
            'Pathway_ID': pw_id,
            'Description': pw_name,
            'GeneRatio': f"{k}/{n}",
            'BgRatio': f"{K}/{N}",
            'pvalue': pvalue,
            'OddsRatio': round(odds_ratio, 4),
            'Count': k,
            'Genes': ';'.join(pw_genes.intersection(de_mapped))
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    # Benjamini-Hochberg correction
    _, padj, _, _ = multipletests(df['pvalue'].values, method='fdr_bh')
    df['p.adjust'] = padj
    df = df.sort_values('pvalue').reset_index(drop=True)
    return df


def run_kegg_gsea(gene_ranking, gene_pathway_links, pathways,
                  id_mapping=None, min_size=5, max_size=500, permutations=1000):
    """
    Perform GSEA for KEGG pathways using prerank approach.

    Parameters
    ----------
    gene_ranking : pd.Series
        Gene ranking (index=gene_id, values=log2FC), sorted descending.
    gene_pathway_links : dict
        {gene_id: [pathway_ids]}
    pathways : dict
        {pathway_id: pathway_name}
    id_mapping : dict, optional
        Mapping from input gene IDs to KEGG gene IDs.
    min_size : int
        Minimum gene set size.
    max_size : int
        Maximum gene set size.
    permutations : int
        Number of permutations.

    Returns
    -------
    pd.DataFrame
        GSEA results.
    """
    import gseapy as gp

    # Map gene IDs
    if id_mapping:
        reverse_map = {v: k for k, v in id_mapping.items()}
        ranking_mapped = gene_ranking.copy()
        new_index = []
        for idx in ranking_mapped.index:
            if idx in id_mapping:
                new_index.append(id_mapping[idx])
            else:
                new_index.append(idx)
        ranking_mapped.index = new_index
    else:
        ranking_mapped = gene_ranking

    # Remove duplicates in ranking
    ranking_mapped = ranking_mapped[~ranking_mapped.index.duplicated(keep='first')]

    # Build pathway gene sets using pathway names
    pathway_gene_sets = {}
    for gene, pw_list in gene_pathway_links.items():
        for pw in pw_list:
            pw_name = pathways.get(pw, pw)
            if pw_name not in pathway_gene_sets:
                pathway_gene_sets[pw_name] = []
            pathway_gene_sets[pw_name].append(gene)

    # Filter by size
    pathway_gene_sets = {
        k: list(set(v)) for k, v in pathway_gene_sets.items()
        if min_size <= len(set(v)) <= max_size
    }

    if not pathway_gene_sets or len(ranking_mapped) < 10:
        return pd.DataFrame()

    try:
        pre_res = gp.prerank(
            rnk=ranking_mapped,
            gene_sets=pathway_gene_sets,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutations,
            outdir=None,
            seed=42,
            verbose=False,
            no_plot=True
        )
        res_df = pre_res.res2d.copy()
        res_df = res_df.sort_values('NOM p-val' if 'NOM p-val' in res_df.columns else 'pval')
        return res_df
    except Exception as e:
        st.warning(f"GSEA computation warning: {str(e)}")
        return pd.DataFrame()


def get_pathway_image_url(pathway_id, org_code):
    """Get the URL for a KEGG pathway map image."""
    return f"https://www.kegg.jp/kegg-bin/show_pathway?{pathway_id}"
