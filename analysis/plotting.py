"""
Plotting module for Functional Enrichment & GSEA App.
Creates publication-quality interactive and static plots.
"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pandas as pd
import numpy as np
import io
import base64

# ── Color Palettes ──────────────────────────────────────────────────────
PALETTE_UP = '#e74c3c'
PALETTE_DOWN = '#3498db'
PALETTE_NS = '#bdc3c7'
PALETTE_ACCENT = '#2ecc71'
GRADIENT_LOW = '#fee0d2'
GRADIENT_HIGH = '#de2d26'
GRADIENT_GSEA_POS = '#fc8d59'
GRADIENT_GSEA_NEG = '#91bfdb'

PLOTLY_TEMPLATE = "plotly_white"

COG_COLORS = {
    "A": "#e6194b", "B": "#3cb44b", "C": "#ffe119", "D": "#4363d8",
    "E": "#f58231", "F": "#911eb4", "G": "#46f0f0", "H": "#f032e6",
    "I": "#bcf60c", "J": "#fabebe", "K": "#008080", "L": "#e6beff",
    "M": "#9a6324", "N": "#fffac8", "O": "#800000", "P": "#aaffc3",
    "Q": "#808000", "R": "#ffd8b1", "S": "#000075", "T": "#808080",
    "U": "#2f4f4f", "V": "#006400", "W": "#d7f9fa", "X": "#8c7088",
    "Y": "#c0c0c0", "Z": "#ff6347",
}


# ═══════════════════════════════════════════════════════════════════════
#  DATA OVERVIEW PLOTS
# ═══════════════════════════════════════════════════════════════════════

def plot_volcano(df, padj_cutoff=0.05, log2fc_cutoff=1.0):
    """Create an interactive volcano plot."""
    plot_df = df.dropna(subset=['log2FoldChange', 'padj']).copy()
    plot_df['-log10(padj)'] = -np.log10(plot_df['padj'].clip(lower=1e-300))

    # Classify genes
    conditions = [
        (plot_df['padj'] < padj_cutoff) & (plot_df['log2FoldChange'] > log2fc_cutoff),
        (plot_df['padj'] < padj_cutoff) & (plot_df['log2FoldChange'] < -log2fc_cutoff),
    ]
    choices = ['Upregulated', 'Downregulated']
    plot_df['Regulation'] = np.select(conditions, choices, default='Not Significant')

    color_map = {
        'Upregulated': PALETTE_UP,
        'Downregulated': PALETTE_DOWN,
        'Not Significant': PALETTE_NS
    }

    hover_cols = ['gene_id']
    if 'gene_name' in plot_df.columns:
        hover_cols.append('gene_name')
    if 'product' in plot_df.columns:
        hover_cols.append('product')

    fig = px.scatter(
        plot_df, x='log2FoldChange', y='-log10(padj)',
        color='Regulation', color_discrete_map=color_map,
        hover_data=hover_cols,
        template=PLOTLY_TEMPLATE,
        labels={'log2FoldChange': 'log₂(Fold Change)', '-log10(padj)': '-log₁₀(adjusted p-value)'},
    )

    # Add threshold lines
    fig.add_hline(y=-np.log10(padj_cutoff), line_dash="dash",
                  line_color="gray", opacity=0.6,
                  annotation_text=f"padj = {padj_cutoff}")
    fig.add_vline(x=log2fc_cutoff, line_dash="dash", line_color="gray", opacity=0.6)
    fig.add_vline(x=-log2fc_cutoff, line_dash="dash", line_color="gray", opacity=0.6)

    n_up = (plot_df['Regulation'] == 'Upregulated').sum()
    n_down = (plot_df['Regulation'] == 'Downregulated').sum()

    fig.update_layout(
        title=dict(text=f"Volcano Plot — {n_up} Up, {n_down} Down", x=0.5),
        font=dict(family="Arial", size=13),
        legend=dict(title="", orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        height=600, width=850,
        margin=dict(t=80, b=60),
    )
    return fig


def plot_ma(df, padj_cutoff=0.05, log2fc_cutoff=1.0):
    """Create an MA plot (log2FC vs baseMean)."""
    plot_df = df.dropna(subset=['log2FoldChange', 'padj', 'baseMean']).copy()
    plot_df['log10_baseMean'] = np.log10(plot_df['baseMean'].clip(lower=0.01))

    conditions = [
        (plot_df['padj'] < padj_cutoff) & (plot_df['log2FoldChange'] > log2fc_cutoff),
        (plot_df['padj'] < padj_cutoff) & (plot_df['log2FoldChange'] < -log2fc_cutoff),
    ]
    choices = ['Upregulated', 'Downregulated']
    plot_df['Regulation'] = np.select(conditions, choices, default='Not Significant')

    color_map = {
        'Upregulated': PALETTE_UP,
        'Downregulated': PALETTE_DOWN,
        'Not Significant': PALETTE_NS
    }

    fig = px.scatter(
        plot_df, x='log10_baseMean', y='log2FoldChange',
        color='Regulation', color_discrete_map=color_map,
        template=PLOTLY_TEMPLATE,
        labels={'log10_baseMean': 'log₁₀(Base Mean)', 'log2FoldChange': 'log₂(Fold Change)'},
    )
    fig.add_hline(y=0, line_color="gray", opacity=0.5)
    fig.update_layout(
        title=dict(text="MA Plot", x=0.5),
        font=dict(family="Arial", size=13),
        legend=dict(title="", orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        height=550, width=850,
    )
    return fig


def plot_pvalue_histogram(df):
    """Plot p-value distribution."""
    plot_df = df.dropna(subset=['pvalue']).copy()
    fig = px.histogram(
        plot_df, x='pvalue', nbins=50,
        template=PLOTLY_TEMPLATE,
        labels={'pvalue': 'p-value', 'count': 'Frequency'},
        color_discrete_sequence=['#3498db'],
    )
    fig.update_layout(
        title=dict(text="P-value Distribution", x=0.5),
        font=dict(family="Arial", size=13),
        height=400, width=700,
        bargap=0.02,
    )
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  KEGG ENRICHMENT PLOTS
# ═══════════════════════════════════════════════════════════════════════

def plot_kegg_dotplot(ora_df, top_n=20, title="KEGG Pathway Enrichment — Dot Plot"):
    """Create a dot plot for KEGG ORA results."""
    if ora_df.empty:
        return _empty_figure("No significant KEGG pathways found")

    df = ora_df.head(top_n).copy()
    # Parse GeneRatio
    df['GeneRatio_val'] = df['GeneRatio'].apply(lambda x: eval(x) if isinstance(x, str) else x)

    fig = px.scatter(
        df,
        x='GeneRatio_val',
        y='Description',
        size='Count',
        color='p.adjust',
        color_continuous_scale='RdYlBu_r',
        template=PLOTLY_TEMPLATE,
        labels={
            'GeneRatio_val': 'Gene Ratio',
            'Description': '',
            'p.adjust': 'Adjusted p-value',
            'Count': 'Gene Count'
        },
        size_max=22,
    )
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=12),
        yaxis=dict(autorange="reversed", tickfont=dict(size=11)),
        coloraxis_colorbar=dict(title="Adj. p-value"),
        height=max(400, top_n * 28 + 100),
        width=900,
        margin=dict(l=300),
    )
    return fig


def plot_kegg_barplot(ora_df, top_n=20, title="KEGG Pathway Enrichment — Bar Plot"):
    """Create a bar plot for KEGG ORA results."""
    if ora_df.empty:
        return _empty_figure("No significant KEGG pathways found")

    df = ora_df.head(top_n).copy()
    df['-log10(p.adjust)'] = -np.log10(df['p.adjust'].clip(lower=1e-300))

    fig = px.bar(
        df,
        x='-log10(p.adjust)',
        y='Description',
        color='Count',
        color_continuous_scale='Viridis',
        template=PLOTLY_TEMPLATE,
        orientation='h',
        labels={
            '-log10(p.adjust)': '-log₁₀(adjusted p-value)',
            'Description': '',
            'Count': 'Gene Count'
        },
    )
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=12),
        yaxis=dict(autorange="reversed", tickfont=dict(size=11)),
        height=max(400, top_n * 28 + 100),
        width=900,
        margin=dict(l=300),
    )
    return fig


def plot_kegg_lollipop(ora_df, top_n=15, title="KEGG Pathway Enrichment — Lollipop Plot"):
    """Create a lollipop plot for KEGG ORA results."""
    if ora_df.empty:
        return _empty_figure("No significant KEGG pathways found")

    df = ora_df.head(top_n).copy()
    df['-log10(p.adjust)'] = -np.log10(df['p.adjust'].clip(lower=1e-300))
    df = df.sort_values('-log10(p.adjust)')

    fig = go.Figure()
    for i, row in df.iterrows():
        fig.add_trace(go.Scatter(
            x=[0, row['-log10(p.adjust)']],
            y=[row['Description'], row['Description']],
            mode='lines',
            line=dict(color='#d3d3d3', width=2),
            showlegend=False,
        ))
    fig.add_trace(go.Scatter(
        x=df['-log10(p.adjust)'],
        y=df['Description'],
        mode='markers',
        marker=dict(
            size=df['Count'] * 2 + 6,
            color=df['p.adjust'],
            colorscale='RdYlBu_r',
            colorbar=dict(title="Adj. p-value"),
            line=dict(width=1, color='#333'),
        ),
        text=df['Count'],
        hovertemplate="<b>%{y}</b><br>-log10(padj): %{x:.2f}<br>Count: %{text}<extra></extra>",
        showlegend=False,
    ))
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title="-log₁₀(adjusted p-value)",
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=12),
        height=max(400, top_n * 30 + 100),
        width=900,
        margin=dict(l=320),
    )
    return fig


def plot_kegg_network(ora_df, top_n=10, title="KEGG Pathway Enrichment — Network"):
    """Create a simple network-style treemap for enriched pathways."""
    if ora_df.empty:
        return _empty_figure("No significant KEGG pathways found")

    df = ora_df.head(top_n).copy()
    fig = px.treemap(
        df,
        path=['Description'],
        values='Count',
        color='p.adjust',
        color_continuous_scale='RdYlBu_r',
        template=PLOTLY_TEMPLATE,
    )
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=13),
        height=500, width=900,
    )
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  GSEA PLOTS
# ═══════════════════════════════════════════════════════════════════════

def plot_gsea_dotplot(gsea_df, top_n=20, title="GSEA — Enrichment Dot Plot"):
    """Dot plot for GSEA results with positive and negative enrichment."""
    if gsea_df.empty:
        return _empty_figure("No significant GSEA results found")

    # Normalize column names
    col_map = {}
    for c in gsea_df.columns:
        cl = c.lower()
        # Check lead_genes/genes BEFORE nes (since 'lead_genes' contains 'nes')
        if 'lead_genes' in cl or 'ledge_genes' in cl:
            col_map[c] = 'Genes'
        elif cl == 'nes' or cl == 'normalized enrichment score':
            col_map[c] = 'NES'
        elif 'nom p-val' in cl or cl == 'pval' or 'nom_p_val' in cl:
            col_map[c] = 'pval'
        elif 'fdr' in cl or 'fdr q-val' in cl:
            col_map[c] = 'fdr'
        elif cl == 'term':
            col_map[c] = 'Term'
        elif 'geneset_size' in cl or cl == 'size':
            col_map[c] = 'Size'
        elif cl == 'es':
            col_map[c] = 'ES'
        elif 'genes' in cl and 'lead' not in cl:
            col_map[c] = 'Genes'

    df = gsea_df.rename(columns=col_map).copy()

    if 'Term' not in df.columns:
        df['Term'] = df.index

    if 'NES' not in df.columns:
        return _empty_figure("GSEA results missing NES column")

    df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
    if 'pval' in df.columns:
        df['pval'] = pd.to_numeric(df['pval'], errors='coerce')
    if 'fdr' in df.columns:
        df['fdr'] = pd.to_numeric(df['fdr'], errors='coerce')

    # Split into positive and negative
    df_pos = df[df['NES'] > 0].nlargest(top_n // 2, 'NES')
    df_neg = df[df['NES'] < 0].nsmallest(top_n // 2, 'NES')
    df_combined = pd.concat([df_pos, df_neg]).sort_values('NES')

    if df_combined.empty:
        return _empty_figure("No significant GSEA results found")

    color_col = 'fdr' if 'fdr' in df_combined.columns else 'pval'
    size_col = 'Size' if 'Size' in df_combined.columns else None

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df_combined['NES'],
        y=df_combined['Term'],
        mode='markers',
        marker=dict(
            size=df_combined[size_col].clip(upper=100) / 3 + 6 if size_col else 14,
            color=df_combined[color_col] if color_col in df_combined.columns else '#3498db',
            colorscale='RdYlBu_r',
            colorbar=dict(title="FDR q-value" if color_col == 'fdr' else "p-value"),
            line=dict(width=1, color='#555'),
        ),
        hovertemplate="<b>%{y}</b><br>NES: %{x:.3f}<extra></extra>",
    ))

    fig.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.5)

    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=12),
        height=max(450, len(df_combined) * 28 + 100),
        width=950,
        margin=dict(l=350),
        showlegend=False,
    )
    return fig


def plot_gsea_barplot(gsea_df, top_n=20, title="GSEA — NES Bar Plot"):
    """Bar plot showing NES values colored by direction."""
    if gsea_df.empty:
        return _empty_figure("No significant GSEA results found")

    col_map = {}
    for c in gsea_df.columns:
        cl = c.lower()
        if 'lead_genes' in cl or 'ledge_genes' in cl:
            continue  # skip gene lists
        elif cl == 'nes':
            col_map[c] = 'NES'
        elif cl == 'term':
            col_map[c] = 'Term'
        elif 'fdr' in cl:
            col_map[c] = 'fdr'
        elif 'nom p-val' in cl or cl == 'pval':
            col_map[c] = 'pval'

    df = gsea_df.rename(columns=col_map).copy()
    if 'Term' not in df.columns:
        df['Term'] = df.index
    if 'NES' not in df.columns:
        return _empty_figure("GSEA results missing NES column")

    df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
    df_pos = df[df['NES'] > 0].nlargest(top_n // 2, 'NES')
    df_neg = df[df['NES'] < 0].nsmallest(top_n // 2, 'NES')
    df_combined = pd.concat([df_pos, df_neg]).sort_values('NES')

    if df_combined.empty:
        return _empty_figure("No significant GSEA results found")

    colors = [PALETTE_UP if v > 0 else PALETTE_DOWN for v in df_combined['NES']]

    fig = go.Figure(go.Bar(
        x=df_combined['NES'],
        y=df_combined['Term'],
        orientation='h',
        marker_color=colors,
        hovertemplate="<b>%{y}</b><br>NES: %{x:.3f}<extra></extra>",
    ))

    fig.add_vline(x=0, line_color="gray", opacity=0.5)

    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=12),
        height=max(400, len(df_combined) * 28 + 100),
        width=900,
        margin=dict(l=350),
    )
    return fig


def plot_gsea_waterfall(gsea_df, top_n=30, title="GSEA — Waterfall Plot"):
    """Waterfall plot showing ranked NES values."""
    if gsea_df.empty:
        return _empty_figure("No significant GSEA results found")

    col_map = {}
    for c in gsea_df.columns:
        cl = c.lower()
        if 'lead_genes' in cl or 'ledge_genes' in cl:
            continue  # skip gene lists
        elif cl == 'nes':
            col_map[c] = 'NES'
        elif cl == 'term':
            col_map[c] = 'Term'

    df = gsea_df.rename(columns=col_map).copy()
    if 'Term' not in df.columns:
        df['Term'] = df.index
    if 'NES' not in df.columns:
        return _empty_figure("GSEA results missing NES column")

    df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
    df = df.dropna(subset=['NES']).sort_values('NES', ascending=False).head(top_n)

    colors = [PALETTE_UP if v > 0 else PALETTE_DOWN for v in df['NES']]

    fig = go.Figure(go.Bar(
        x=list(range(len(df))),
        y=df['NES'],
        marker_color=colors,
        text=df['Term'],
        hovertemplate="<b>%{text}</b><br>NES: %{y:.3f}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(text=title, x=0.5),
        yaxis_title="Normalized Enrichment Score (NES)",
        xaxis_title="Ranked Pathways",
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=12),
        height=500, width=900,
        xaxis=dict(showticklabels=False),
    )
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  COG ANALYSIS PLOTS
# ═══════════════════════════════════════════════════════════════════════

def plot_cog_distribution_bar(dist_dfs, title="COG Category Distribution"):
    """
    Grouped bar chart showing COG category distribution for multiple groups.

    Parameters
    ----------
    dist_dfs : list of pd.DataFrame
        Each with columns: COG_Category, Count, Percentage, Group
    """
    if not dist_dfs:
        return _empty_figure("No COG distribution data")

    combined = pd.concat(dist_dfs, ignore_index=True)
    combined = combined[combined['Count'] > 0]

    if combined.empty:
        return _empty_figure("No COG category assignments found")

    fig = px.bar(
        combined,
        x='COG_Category',
        y='Percentage',
        color='Group',
        barmode='group',
        template=PLOTLY_TEMPLATE,
        color_discrete_sequence=[PALETTE_NS, PALETTE_UP, PALETTE_DOWN],
        labels={'COG_Category': 'COG Category', 'Percentage': 'Percentage (%)'},
        hover_data=['Description', 'Count'],
    )
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=13),
        legend=dict(title="", orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        height=550, width=1000,
        xaxis=dict(dtick=1),
    )
    return fig


def plot_cog_distribution_stacked(dist_dfs, title="COG Category Distribution — Stacked"):
    """Stacked bar chart of COG category percentages per group."""
    if not dist_dfs:
        return _empty_figure("No COG distribution data")

    combined = pd.concat(dist_dfs, ignore_index=True)
    combined = combined[combined['Count'] > 0]
    if combined.empty:
        return _empty_figure("No COG data")

    fig = px.bar(
        combined,
        x='Group',
        y='Percentage',
        color='COG_Category',
        template=PLOTLY_TEMPLATE,
        color_discrete_map=COG_COLORS,
        hover_data=['Description', 'Count'],
        labels={'Percentage': 'Percentage (%)', 'COG_Category': 'COG'},
    )
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=13),
        barmode='stack',
        height=550, width=800,
        legend=dict(title="COG Category"),
    )
    return fig


def plot_cog_pie(dist_df, title="COG Category Proportions"):
    """Pie chart showing COG category proportions."""
    df = dist_df[dist_df['Count'] > 0].copy()
    if df.empty:
        return _empty_figure("No COG data")

    colors = [COG_COLORS.get(c, '#999') for c in df['COG_Category']]

    fig = go.Figure(go.Pie(
        labels=df['COG_Category'] + ': ' + df['Description'],
        values=df['Count'],
        marker=dict(colors=colors),
        textinfo='label+percent',
        textposition='outside',
        hole=0.35,
    ))
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=11),
        height=650, width=900,
        showlegend=True,
        legend=dict(font=dict(size=10)),
    )
    return fig


def plot_cog_enrichment_dotplot(enrich_df, top_n=20,
                                 title="COG Category Enrichment — Dot Plot"):
    """Dot plot for COG enrichment results."""
    if enrich_df.empty:
        return _empty_figure("No enriched COG categories found")

    df = enrich_df[enrich_df['DE_Count'] > 0].head(top_n).copy()
    if df.empty:
        return _empty_figure("No enriched COG categories found")

    df['Label'] = df['COG_Category'] + ': ' + df['Description']
    df['-log10(p.adjust)'] = -np.log10(df['p.adjust'].clip(lower=1e-300))

    fig = px.scatter(
        df,
        x='Fold_Enrichment',
        y='Label',
        size='DE_Count',
        color='p.adjust',
        color_continuous_scale='RdYlBu_r',
        template=PLOTLY_TEMPLATE,
        labels={
            'Fold_Enrichment': 'Fold Enrichment',
            'Label': '',
            'p.adjust': 'Adjusted p-value',
            'DE_Count': 'DE Gene Count'
        },
        size_max=20,
    )
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)
    fig.update_layout(
        title=dict(text=title, x=0.5),
        font=dict(family="Arial", size=12),
        yaxis=dict(autorange="reversed"),
        height=max(400, len(df) * 30 + 100),
        width=900,
        margin=dict(l=350),
    )
    return fig


def plot_cog_enrichment_bar(enrich_df, top_n=20,
                             title="COG Category Enrichment — Bar Plot"):
    """Bar plot for COG enrichment results."""
    if enrich_df.empty:
        return _empty_figure("No enriched COG categories found")

    df = enrich_df[enrich_df['DE_Count'] > 0].head(top_n).copy()
    if df.empty:
        return _empty_figure("No enriched COG categories found")

    df['Label'] = df['COG_Category'] + ': ' + df['Description']
    df['-log10(p.adjust)'] = -np.log10(df['p.adjust'].clip(lower=1e-300))
    colors = [COG_COLORS.get(c, '#999') for c in df['COG_Category']]

    fig = go.Figure(go.Bar(
        x=df['-log10(p.adjust)'],
        y=df['Label'],
        orientation='h',
        marker_color=colors,
        text=df['DE_Count'],
        textposition='outside',
        hovertemplate="<b>%{y}</b><br>-log10(padj): %{x:.2f}<br>Count: %{text}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title="-log₁₀(adjusted p-value)",
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=12),
        yaxis=dict(autorange="reversed"),
        height=max(400, len(df) * 30 + 100),
        width=900,
        margin=dict(l=350),
    )
    return fig


def plot_cog_heatmap(up_dist, down_dist, title="COG Category Heatmap — Up vs Down"):
    """Heatmap showing COG category percentages for up vs down regulated genes."""
    if up_dist.empty and down_dist.empty:
        return _empty_figure("No COG data")

    # Merge
    merged = pd.merge(
        up_dist[['COG_Category', 'Description', 'Percentage']].rename(columns={'Percentage': 'Up (%)'}),
        down_dist[['COG_Category', 'Description', 'Percentage']].rename(columns={'Percentage': 'Down (%)'}),
        on=['COG_Category', 'Description'], how='outer'
    ).fillna(0)

    merged = merged[(merged['Up (%)'] > 0) | (merged['Down (%)'] > 0)]
    if merged.empty:
        return _empty_figure("No COG data")

    merged['Label'] = merged['COG_Category'] + ': ' + merged['Description']

    z_data = merged[['Up (%)', 'Down (%)']].values
    fig = go.Figure(go.Heatmap(
        z=z_data,
        y=merged['Label'],
        x=['Upregulated', 'Downregulated'],
        colorscale='RdBu_r',
        texttemplate='%{z:.1f}%',
        hovertemplate="<b>%{y}</b><br>%{x}: %{z:.1f}%<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=title, x=0.5),
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=12),
        height=max(400, len(merged) * 25 + 100),
        width=700,
        margin=dict(l=350),
    )
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  COMBINED / COMPARISON PLOTS
# ═══════════════════════════════════════════════════════════════════════

def plot_kegg_cog_summary(kegg_ora_df, cog_enrich_df, title="KEGG vs COG Enrichment Summary"):
    """Side-by-side comparison of top KEGG and COG enrichments."""
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Top KEGG Pathways", "Top COG Categories"),
                        horizontal_spacing=0.35)

    # KEGG subplot
    if not kegg_ora_df.empty:
        kegg_top = kegg_ora_df.head(10)
        kegg_top['-log10(padj)'] = -np.log10(kegg_top['p.adjust'].clip(lower=1e-300))
        fig.add_trace(go.Bar(
            x=kegg_top['-log10(padj)'],
            y=kegg_top['Description'],
            orientation='h',
            marker_color='#3498db',
            name='KEGG',
            showlegend=False,
        ), row=1, col=1)

    # COG subplot
    if not cog_enrich_df.empty:
        cog_top = cog_enrich_df[cog_enrich_df['DE_Count'] > 0].head(10)
        if not cog_top.empty:
            cog_top['-log10(padj)'] = -np.log10(cog_top['p.adjust'].clip(lower=1e-300))
            labels = cog_top['COG_Category'] + ': ' + cog_top['Description'].str[:30]
            colors = [COG_COLORS.get(c, '#999') for c in cog_top['COG_Category']]
            fig.add_trace(go.Bar(
                x=cog_top['-log10(padj)'],
                y=labels,
                orientation='h',
                marker_color=colors,
                name='COG',
                showlegend=False,
            ), row=1, col=2)

    fig.update_layout(
        title=dict(text=title, x=0.5),
        template=PLOTLY_TEMPLATE,
        font=dict(family="Arial", size=11),
        height=550, width=1200,
    )
    fig.update_yaxes(autorange="reversed")
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  UTILITY
# ═══════════════════════════════════════════════════════════════════════

def _empty_figure(message="No data to display"):
    """Create a placeholder figure with a message."""
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper", yref="paper",
        x=0.5, y=0.5,
        showarrow=False,
        font=dict(size=16, color="gray"),
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        height=300, width=600,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    return fig


def fig_to_bytes(fig, format='png', scale=3):
    """Convert a Plotly figure to bytes for download."""
    return fig.to_image(format=format, scale=scale)
