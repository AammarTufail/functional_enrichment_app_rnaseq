"""
╔══════════════════════════════════════════════════════════════════════╗
║   Functional Enrichment & GSEA Analysis App                        ║
║   ─────────────────────────────────────────                        ║
║   KEGG Pathway Enrichment · COG Functional Analysis · GSEA        ║
║                                                                    ║
║   Inspired by:  github.com/RickGelhausen/pathsnake                 ║
║   Built with:   Streamlit · Plotly · gseapy · scipy                ║
╚══════════════════════════════════════════════════════════════════════╝
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import time
import warnings
warnings.filterwarnings('ignore')

# ── Page Configuration ────────────────────────────────────────────────
st.set_page_config(
    page_title="Functional Enrichment & GSEA",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────────────
st.markdown("""
<style>
    /* Main header */
    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 2rem 2.5rem;
        border-radius: 12px;
        margin-bottom: 1.5rem;
        color: #ffffff;
        text-align: center;
        border: 1px solid rgba(255,255,255,0.1);
        box-shadow: 0 4px 20px rgba(0,0,0,0.3);
    }
    .main-header h1 {
        font-size: 2.2rem;
        font-weight: 700;
        margin-bottom: 0.3rem;
        letter-spacing: 0.5px;
    }
    .main-header p {
        font-size: 1rem;
        opacity: 0.85;
        margin: 0;
    }

    /* Metric cards */
    .metric-card {
        background: #f8f9fa;
        border-radius: 10px;
        padding: 1.2rem;
        text-align: center;
        border: 1px solid #e9ecef;
        transition: transform 0.2s;
    }
    .metric-card:hover { transform: translateY(-2px); }
    .metric-card h3 {
        font-size: 2rem;
        color: #2c3e50;
        margin: 0;
        font-weight: 700;
    }
    .metric-card p {
        font-size: 0.85rem;
        color: #7f8c8d;
        margin: 0.3rem 0 0 0;
    }
    .metric-up h3 { color: #e74c3c; }
    .metric-down h3 { color: #3498db; }
    .metric-total h3 { color: #2ecc71; }

    /* Section headers */
    .section-header {
        background: linear-gradient(90deg, #f8f9fa 0%, white 100%);
        padding: 0.8rem 1.2rem;
        border-left: 4px solid #3498db;
        border-radius: 0 8px 8px 0;
        margin: 1.5rem 0 1rem 0;
        font-weight: 600;
        font-size: 1.1rem;
        color: #2c3e50;
    }

    /* Info panels */
    .info-panel {
        background: #eef6ff;
        border: 1px solid #bdd7ee;
        border-radius: 8px;
        padding: 1rem 1.2rem;
        margin: 0.8rem 0;
        font-size: 0.9rem;
        color: #1a5276;
    }

    /* Sidebar styling */
    section[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #fafbfc 0%, #f0f2f6 100%);
    }
    /* Sidebar text colors — labels, headings, markdown */
    section[data-testid="stSidebar"] h1,
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3,
    section[data-testid="stSidebar"] h4,
    section[data-testid="stSidebar"] label,
    section[data-testid="stSidebar"] p,
    section[data-testid="stSidebar"] .stMarkdown,
    section[data-testid="stSidebar"] .stMarkdown span,
    section[data-testid="stSidebar"] .stMarkdown div,
    section[data-testid="stSidebar"] .stRadio label,
    section[data-testid="stSidebar"] .stCheckbox label,
    section[data-testid="stSidebar"] .stSelectbox label,
    section[data-testid="stSidebar"] .stSlider label,
    section[data-testid="stSidebar"] [data-testid="stWidgetLabel"],
    section[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {
        color: #1a1a2e !important;
    }
    section[data-testid="stSidebar"] .stRadio > label {
        font-weight: 600;
        color: #1a1a2e !important;
    }
    /* Sidebar file uploader — light background so text is visible */
    section[data-testid="stSidebar"] [data-testid="stFileUploader"],
    section[data-testid="stSidebar"] [data-testid="stFileUploader"] div,
    section[data-testid="stSidebar"] [data-testid="stFileUploader"] section,
    section[data-testid="stSidebar"] [data-testid="stFileUploader"] span,
    section[data-testid="stSidebar"] [data-testid="stFileUploader"] small,
    section[data-testid="stSidebar"] [data-testid="stFileUploader"] p,
    section[data-testid="stSidebar"] [data-testid="stFileUploader"] label {
        color: #1a1a2e !important;
    }
    section[data-testid="stSidebar"] [data-testid="stFileUploaderDropzone"] {
        background-color: #ffffff !important;
        border: 2px dashed #b0b8c4 !important;
        color: #1a1a2e !important;
    }
    section[data-testid="stSidebar"] [data-testid="stFileUploaderDropzone"] * {
        color: #1a1a2e !important;
    }
    /* Sidebar text inputs — force light background */
    section[data-testid="stSidebar"] input[type="text"],
    section[data-testid="stSidebar"] .stTextInput input {
        background-color: #ffffff !important;
        color: #1a1a2e !important;
        border: 1px solid #b0b8c4 !important;
    }
    /* Sidebar select boxes */
    section[data-testid="stSidebar"] [data-baseweb="select"] {
        background-color: #ffffff !important;
    }
    section[data-testid="stSidebar"] [data-baseweb="select"] * {
        color: #1a1a2e !important;
    }
    /* Sidebar slider values and thumb */
    section[data-testid="stSidebar"] [data-testid="stThumbValue"],
    section[data-testid="stSidebar"] [data-testid="stTickBarMin"],
    section[data-testid="stSidebar"] [data-testid="stTickBarMax"] {
        color: #1a1a2e !important;
    }
    /* Sidebar checkbox */
    section[data-testid="stSidebar"] .stCheckbox span {
        color: #1a1a2e !important;
    }
    /* Sidebar links */
    section[data-testid="stSidebar"] .stMarkdown a {
        color: #3498db !important;
    }
    section[data-testid="stSidebar"] hr {
        border-color: #ccc !important;
    }
    /* Sidebar button text should stay white */
    section[data-testid="stSidebar"] button[kind="primary"],
    section[data-testid="stSidebar"] button[kind="primary"] p {
        color: #ffffff !important;
    }

    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px 8px 0 0;
        padding: 10px 24px;
        font-weight: 600;
    }

    /* Download buttons */
    .stDownloadButton button {
        background: linear-gradient(135deg, #2ecc71, #27ae60);
        color: white;
        border: none;
        border-radius: 8px;
        font-weight: 600;
    }

    /* Hide default header decoration */
    header { visibility: hidden; }
    .block-container { padding-top: 1rem; }
</style>
""", unsafe_allow_html=True)


# ── Imports for analysis ──────────────────────────────────────────────
from analysis.data_parser import parse_deseq2_results, get_de_genes, get_gene_ranking
from analysis.kegg_analysis import (
    fetch_kegg_gene_list, fetch_kegg_pathways,
    fetch_kegg_gene_pathway_links, build_kegg_id_mapping,
    run_kegg_ora, run_kegg_gsea, get_pathway_image_url
)
from analysis.cog_analysis import (
    COG_CATEGORIES, COG_SUPERCATEGORIES, COG_COLORS,
    fetch_cog_from_kegg, parse_cog_annotation_file,
    infer_cog_from_products, run_cog_enrichment, get_cog_distribution
)
from analysis.plotting import (
    plot_volcano, plot_ma, plot_pvalue_histogram,
    plot_kegg_dotplot, plot_kegg_barplot, plot_kegg_lollipop, plot_kegg_network,
    plot_gsea_dotplot, plot_gsea_barplot, plot_gsea_waterfall,
    plot_cog_distribution_bar, plot_cog_distribution_stacked,
    plot_cog_pie, plot_cog_enrichment_dotplot, plot_cog_enrichment_bar,
    plot_cog_heatmap, plot_kegg_cog_summary,
)


# ═══════════════════════════════════════════════════════════════════════
#  HEADER
# ═══════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="main-header">
    <h1>🧬 Functional Enrichment & GSEA Analysis</h1>
    <p>KEGG Pathway Enrichment · COG Functional Classification · Gene Set Enrichment Analysis</p>
</div>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════
#  SIDEBAR — Configuration
# ═══════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown("## ⚙️ Configuration")
    st.markdown("---")

    # ── File Upload ───────────────────────────────────────────────────
    st.markdown("### 📁 Data Input")
    uploaded_file = st.file_uploader(
        "Upload DESeq2 Results",
        type=['csv', 'tsv', 'txt'],
        help="Tab or comma-separated file with log2FoldChange and padj columns"
    )

    use_example = st.checkbox("Use example data", value=False,
                               help="Load the included example DESeq2 results")

    st.markdown("---")

    # ── Organism Settings ─────────────────────────────────────────────
    st.markdown("### 🦠 Organism")
    st.markdown(
        "Find your organism code at "
        "[KEGG Organism List](https://www.genome.jp/kegg/catalog/org_list.html)",
        unsafe_allow_html=False,
    )
    kegg_org_code = st.text_input(
        "KEGG Organism Code",
        value="sey",
        help="Three-letter KEGG organism code (e.g., eco, sey, stm). Use 'sey' for SL1344."
    )

    st.markdown("---")

    # ── Analysis Parameters ───────────────────────────────────────────
    st.markdown("### 📊 Parameters")
    padj_cutoff = st.slider("Adjusted p-value cutoff", 0.001, 0.1, 0.05, 0.005,
                             format="%.3f")
    log2fc_cutoff = st.slider("log₂FC cutoff", 0.0, 4.0, 1.0, 0.25)
    top_n_pathways = st.slider("Top pathways to display", 5, 50, 20, 5)
    min_gene_set_size = st.slider("Min gene set size", 2, 15, 3, 1,
                                   help="Minimum genes per pathway for ORA/GSEA")
    gsea_permutations = st.select_slider("GSEA permutations",
                                          options=[100, 500, 1000, 5000, 10000],
                                          value=1000)

    st.markdown("---")

    # ── COG Settings ──────────────────────────────────────────────────
    st.markdown("### 🏷️ COG Analysis")
    cog_source = st.radio(
        "COG annotation source",
        ["Infer from product descriptions", "Fetch from KEGG", "Upload COG file"],
        help="Choose how to obtain COG category assignments"
    )

    cog_file = None
    if cog_source == "Upload COG file":
        cog_file = st.file_uploader(
            "Upload COG Annotation",
            type=['tsv', 'txt', 'csv'],
            help="Tab-separated: gene_id <TAB> COG_categories (e.g., J, KL, E)"
        )

    st.markdown("---")

    # ── Analysis Buttons ──────────────────────────────────────────────
    run_analysis = st.button("🚀 Run Full Analysis", type="primary", use_container_width=True)
    st.markdown("")
    st.markdown(
        "<div style='text-align:center; font-size:0.75rem; color:#999; margin-top:1rem;'>"
        "Inspired by <a href='https://github.com/RickGelhausen/pathsnake' target='_blank'>pathsnake</a>"
        "</div>",
        unsafe_allow_html=True,
    )


# ═══════════════════════════════════════════════════════════════════════
#  LOAD DATA
# ═══════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def load_data(file_input, is_path=False):
    """Load and parse the input file."""
    if is_path:
        return parse_deseq2_results(file_input, sep='\t')
    else:
        return parse_deseq2_results(file_input)


data_loaded = False
df = None

if uploaded_file is not None:
    with st.spinner("Parsing uploaded file..."):
        df = load_data(uploaded_file)
    data_loaded = True
elif use_example:
    import os
    example_path = os.path.join(os.path.dirname(__file__),
                                 "deseq_comp_InSPI2_vs_LSP_with_annotation_and_countings.csv")
    if os.path.exists(example_path):
        with st.spinner("Loading example data..."):
            df = load_data(example_path, is_path=True)
        data_loaded = True
    else:
        st.error("Example file not found. Please upload a file.")

if not data_loaded:
    # Landing page
    st.markdown("""
    <div class="info-panel" style="text-align:center; padding:2rem;">
        <h3 style="margin-top:0;">👋 Welcome!</h3>
        <p>Upload your <b>DESeq2 results file</b> or check <b>"Use example data"</b> to get started.</p>
        <p style="margin-top:1rem;">This app performs:</p>
        <ul style="text-align:left; max-width:400px; margin:0.5rem auto;">
            <li><b>KEGG Pathway Enrichment</b> — Over-Representation Analysis</li>
            <li><b>KEGG GSEA</b> — Gene Set Enrichment Analysis (prerank)</li>
            <li><b>COG Functional Classification</b> — Category enrichment & distribution</li>
        </ul>
        <p style="margin-top:1rem; font-size:0.85rem; color:#999;">
            Supported input: Tab/comma-separated files with log2FoldChange and padj columns
        </p>
    </div>
    """, unsafe_allow_html=True)
    st.stop()


# ═══════════════════════════════════════════════════════════════════════
#  DATA OVERVIEW TAB SYSTEM
# ═══════════════════════════════════════════════════════════════════════

# Get DE genes
up_genes, down_genes, all_de_genes = get_de_genes(df, padj_cutoff, log2fc_cutoff)

# Create tabs
tab_overview, tab_kegg, tab_gsea, tab_cog, tab_combined = st.tabs([
    "📋 Data Overview",
    "🗺️ KEGG Enrichment",
    "📈 KEGG GSEA",
    "🏷️ COG Analysis",
    "📊 Summary",
])


# ═══════════════════════════════════════════════════════════════════════
#  TAB 1: DATA OVERVIEW
# ═══════════════════════════════════════════════════════════════════════
with tab_overview:
    st.markdown('<div class="section-header">📋 Dataset Summary</div>', unsafe_allow_html=True)

    # Metric cards
    total_genes = len(df)
    valid_genes = df.dropna(subset=['log2FoldChange', 'padj']).shape[0]
    n_up = len(up_genes)
    n_down = len(down_genes)

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown(f"""<div class="metric-card metric-total">
            <h3>{total_genes:,}</h3><p>Total Genes</p></div>""", unsafe_allow_html=True)
    with col2:
        st.markdown(f"""<div class="metric-card">
            <h3>{valid_genes:,}</h3><p>Tested Genes</p></div>""", unsafe_allow_html=True)
    with col3:
        st.markdown(f"""<div class="metric-card metric-up">
            <h3>{n_up:,}</h3><p>Upregulated</p></div>""", unsafe_allow_html=True)
    with col4:
        st.markdown(f"""<div class="metric-card metric-down">
            <h3>{n_down:,}</h3><p>Downregulated</p></div>""", unsafe_allow_html=True)

    st.markdown("")

    # Plots
    col_v, col_m = st.columns(2)
    with col_v:
        fig_volcano = plot_volcano(df, padj_cutoff, log2fc_cutoff)
        st.plotly_chart(fig_volcano, use_container_width=True)
    with col_m:
        fig_ma = plot_ma(df, padj_cutoff, log2fc_cutoff)
        st.plotly_chart(fig_ma, use_container_width=True)

    # P-value histogram
    if 'pvalue' in df.columns:
        fig_phist = plot_pvalue_histogram(df)
        st.plotly_chart(fig_phist, use_container_width=True)

    # Data table
    with st.expander("📄 View Data Table", expanded=False):
        display_cols = ['gene_id']
        if 'gene_name' in df.columns:
            display_cols.append('gene_name')
        if 'product' in df.columns:
            display_cols.append('product')
        display_cols.extend(['log2FoldChange', 'padj', 'pvalue', 'baseMean'])
        display_cols = [c for c in display_cols if c in df.columns]
        st.dataframe(
            df[display_cols].sort_values('padj').head(500),
            use_container_width=True,
            height=400,
        )

    # Download processed data
    csv_buf = io.BytesIO()
    df.to_csv(csv_buf, index=False, sep='\t')
    st.download_button(
        "⬇️ Download Processed Data (TSV)",
        csv_buf.getvalue(),
        "processed_deseq2_results.tsv",
        "text/tab-separated-values",
    )


# ═══════════════════════════════════════════════════════════════════════
#  ANALYSIS STATE
# ═══════════════════════════════════════════════════════════════════════
# Use session state for analysis results
if 'kegg_ora_up' not in st.session_state:
    st.session_state.kegg_ora_up = pd.DataFrame()
if 'kegg_ora_down' not in st.session_state:
    st.session_state.kegg_ora_down = pd.DataFrame()
if 'kegg_ora_all' not in st.session_state:
    st.session_state.kegg_ora_all = pd.DataFrame()
if 'kegg_gsea_results' not in st.session_state:
    st.session_state.kegg_gsea_results = pd.DataFrame()
if 'cog_mapping' not in st.session_state:
    st.session_state.cog_mapping = {}
if 'cog_enrich_up' not in st.session_state:
    st.session_state.cog_enrich_up = pd.DataFrame()
if 'cog_enrich_down' not in st.session_state:
    st.session_state.cog_enrich_down = pd.DataFrame()
if 'cog_enrich_all' not in st.session_state:
    st.session_state.cog_enrich_all = pd.DataFrame()
if 'id_mapping' not in st.session_state:
    st.session_state.id_mapping = {}
if 'analysis_done' not in st.session_state:
    st.session_state.analysis_done = False


# ═══════════════════════════════════════════════════════════════════════
#  RUN ANALYSIS
# ═══════════════════════════════════════════════════════════════════════
if run_analysis and data_loaded:
    all_gene_ids = df['gene_id'].dropna().unique().tolist()
    up_ids = up_genes['gene_id'].dropna().tolist()
    down_ids = down_genes['gene_id'].dropna().tolist()
    all_de_ids = all_de_genes['gene_id'].dropna().tolist()

    progress = st.progress(0, text="Starting analysis...")

    # ── Step 1: Fetch KEGG Data ────────────────────────────────────────
    progress.progress(5, text="Fetching KEGG gene list...")
    try:
        kegg_genes = fetch_kegg_gene_list(kegg_org_code)
        progress.progress(15, text=f"Found {len(kegg_genes)} KEGG genes. Fetching pathways...")

        kegg_pathways = fetch_kegg_pathways(kegg_org_code)
        progress.progress(25, text=f"Found {len(kegg_pathways)} pathways. Fetching links...")

        gene_pathway_links = fetch_kegg_gene_pathway_links(kegg_org_code)
        progress.progress(35, text=f"Got {len(gene_pathway_links)} gene-pathway links. Building ID mapping...")

        # ── Step 2: Map Gene IDs ───────────────────────────────────────
        # Build gene_name and product maps for cross-strain ID mapping
        gene_name_map = {}
        product_map = {}
        for _, row in df.iterrows():
            gid = row.get('gene_id')
            if gid and not pd.isna(gid):
                gname = row.get('gene_name')
                if gname and not pd.isna(gname):
                    gene_name_map[str(gid)] = str(gname)
                prod = row.get('product')
                if prod and not pd.isna(prod):
                    product_map[str(gid)] = str(prod)
        id_mapping = build_kegg_id_mapping(
            all_gene_ids, kegg_genes, kegg_org_code,
            gene_name_map=gene_name_map, product_map=product_map
        )
        st.session_state.id_mapping = id_mapping
        mapped_pct = len(id_mapping) / len(all_gene_ids) * 100 if all_gene_ids else 0
        progress.progress(40, text=f"Mapped {len(id_mapping)}/{len(all_gene_ids)} genes ({mapped_pct:.1f}%). Running ORA...")

        # ── Step 3: KEGG ORA ───────────────────────────────────────────
        st.session_state.kegg_ora_up = run_kegg_ora(
            up_ids, all_gene_ids, gene_pathway_links, kegg_pathways,
            id_mapping, min_size=min_gene_set_size
        )
        progress.progress(50, text="ORA (upregulated) done. Running ORA (downregulated)...")

        st.session_state.kegg_ora_down = run_kegg_ora(
            down_ids, all_gene_ids, gene_pathway_links, kegg_pathways,
            id_mapping, min_size=min_gene_set_size
        )
        progress.progress(55, text="ORA (downregulated) done. Running ORA (all DE)...")

        st.session_state.kegg_ora_all = run_kegg_ora(
            all_de_ids, all_gene_ids, gene_pathway_links, kegg_pathways,
            id_mapping, min_size=min_gene_set_size
        )
        progress.progress(60, text="ORA done. Running GSEA...")

        # ── Step 4: KEGG GSEA ──────────────────────────────────────────
        gene_ranking = get_gene_ranking(df)
        st.session_state.kegg_gsea_results = run_kegg_gsea(
            gene_ranking, gene_pathway_links, kegg_pathways,
            id_mapping, min_size=min_gene_set_size,
            permutations=gsea_permutations
        )
        progress.progress(75, text="GSEA done. Running COG analysis...")

    except Exception as e:
        st.error(f"KEGG analysis error: {str(e)}")
        progress.progress(75, text="KEGG analysis had issues. Continuing with COG...")

    # ── Step 5: COG Analysis ───────────────────────────────────────────
    try:
        if cog_source == "Upload COG file" and cog_file is not None:
            cog_mapping = parse_cog_annotation_file(cog_file)
        elif cog_source == "Fetch from KEGG":
            progress.progress(80, text="Fetching COG data from KEGG...")
            cog_mapping = fetch_cog_from_kegg(kegg_org_code)
            if not cog_mapping:
                st.info("No COG data from KEGG. Falling back to product-based inference.")
                cog_mapping = infer_cog_from_products(df)
        else:
            cog_mapping = infer_cog_from_products(df)

        st.session_state.cog_mapping = cog_mapping
        progress.progress(85, text=f"COG: {len(cog_mapping)} genes mapped. Running enrichment...")

        # COG Enrichment
        st.session_state.cog_enrich_up = run_cog_enrichment(
            up_ids, all_gene_ids, cog_mapping
        )
        st.session_state.cog_enrich_down = run_cog_enrichment(
            down_ids, all_gene_ids, cog_mapping
        )
        st.session_state.cog_enrich_all = run_cog_enrichment(
            all_de_ids, all_gene_ids, cog_mapping
        )
        progress.progress(95, text="COG enrichment done. Generating plots...")

    except Exception as e:
        st.error(f"COG analysis error: {str(e)}")

    st.session_state.analysis_done = True
    progress.progress(100, text="✅ Analysis complete!")
    time.sleep(0.5)
    progress.empty()
    st.success("Analysis complete! Navigate through the tabs to explore results.")
    st.balloons()


# ═══════════════════════════════════════════════════════════════════════
#  TAB 2: KEGG ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════
with tab_kegg:
    st.markdown('<div class="section-header">🗺️ KEGG Pathway Over-Representation Analysis</div>',
                unsafe_allow_html=True)

    if not st.session_state.analysis_done:
        st.info("👈 Click **Run Full Analysis** in the sidebar to start.")
    else:
        # ID Mapping info
        mapping = st.session_state.id_mapping
        n_mapped = len(mapping)
        n_total = len(df['gene_id'].dropna().unique())
        if n_mapped > 0:
            st.markdown(f"""
            <div class="info-panel">
                <b>Gene ID Mapping:</b> {n_mapped}/{n_total} genes mapped to KEGG
                ({n_mapped/n_total*100:.1f}%) — Organism: <code>{kegg_org_code}</code>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.warning(
                f"Could not map any gene IDs to KEGG organism '{kegg_org_code}'. "
                "Check that the organism code is correct and gene IDs match KEGG entries."
            )

        # Direction selector
        kegg_direction = st.radio(
            "Select gene set",
            ["All DE genes", "Upregulated only", "Downregulated only"],
            horizontal=True, key="kegg_dir"
        )

        if kegg_direction == "All DE genes":
            kegg_ora_df = st.session_state.kegg_ora_all
        elif kegg_direction == "Upregulated only":
            kegg_ora_df = st.session_state.kegg_ora_up
        else:
            kegg_ora_df = st.session_state.kegg_ora_down

        if kegg_ora_df.empty:
            st.warning("No enriched KEGG pathways found for this selection.")
        else:
            # Filter by adjusted p-value
            sig_kegg = kegg_ora_df[kegg_ora_df['p.adjust'] < padj_cutoff]
            st.markdown(f"**{len(sig_kegg)} significant pathways** (padj < {padj_cutoff}), "
                        f"**{len(kegg_ora_df)} total tested**")

            # Plots
            kegg_plot_type = st.selectbox(
                "Plot type",
                ["Dot Plot", "Bar Plot", "Lollipop Plot", "Treemap"],
                key="kegg_plot_type"
            )

            direction_label = kegg_direction.replace(" only", "").replace(" genes", "")
            if kegg_plot_type == "Dot Plot":
                fig = plot_kegg_dotplot(kegg_ora_df, top_n_pathways,
                                        f"KEGG ORA — {direction_label}")
            elif kegg_plot_type == "Bar Plot":
                fig = plot_kegg_barplot(kegg_ora_df, top_n_pathways,
                                        f"KEGG ORA — {direction_label}")
            elif kegg_plot_type == "Lollipop Plot":
                fig = plot_kegg_lollipop(kegg_ora_df, top_n_pathways,
                                          f"KEGG ORA — {direction_label}")
            else:
                fig = plot_kegg_network(kegg_ora_df, top_n_pathways,
                                        f"KEGG ORA — {direction_label}")

            st.plotly_chart(fig, use_container_width=True)

            # KEGG results table
            with st.expander("📄 KEGG ORA Results Table", expanded=False):
                st.dataframe(kegg_ora_df, use_container_width=True, height=400)

            # Download
            col_dl1, col_dl2 = st.columns(2)
            with col_dl1:
                csv_ora = kegg_ora_df.to_csv(index=False, sep='\t')
                st.download_button(
                    "⬇️ Download KEGG ORA Results (TSV)",
                    csv_ora, f"kegg_ora_{kegg_direction.lower().replace(' ', '_')}.tsv",
                    "text/tab-separated-values", key="dl_kegg_ora"
                )
            with col_dl2:
                # Significant only
                if not sig_kegg.empty:
                    csv_sig = sig_kegg.to_csv(index=False, sep='\t')
                    st.download_button(
                        "⬇️ Download Significant Only (TSV)",
                        csv_sig, f"kegg_ora_significant_{kegg_direction.lower().replace(' ', '_')}.tsv",
                        "text/tab-separated-values", key="dl_kegg_sig"
                    )

        # Pathway links
        if not kegg_ora_df.empty:
            with st.expander("🔗 KEGG Pathway Links", expanded=False):
                for _, row in kegg_ora_df.head(10).iterrows():
                    url = get_pathway_image_url(row['Pathway_ID'], kegg_org_code)
                    st.markdown(f"- [{row['Description']}]({url}) (padj={row['p.adjust']:.4e})")


# ═══════════════════════════════════════════════════════════════════════
#  TAB 3: GSEA
# ═══════════════════════════════════════════════════════════════════════
with tab_gsea:
    st.markdown('<div class="section-header">📈 Gene Set Enrichment Analysis (GSEA)</div>',
                unsafe_allow_html=True)

    if not st.session_state.analysis_done:
        st.info("👈 Click **Run Full Analysis** in the sidebar to start.")
    else:
        gsea_df = st.session_state.kegg_gsea_results

        if gsea_df.empty:
            st.warning("No GSEA results available. This may happen when:\n"
                       "- Too few genes mapped to KEGG pathways\n"
                       "- Gene set sizes outside the min/max range\n"
                       "- No significant enrichment detected")
        else:
            st.markdown(f"**{len(gsea_df)} pathways analyzed** via GSEA prerank")

            gsea_plot_type = st.selectbox(
                "Plot type",
                ["Dot Plot", "NES Bar Plot", "Waterfall Plot"],
                key="gsea_plot_type"
            )

            if gsea_plot_type == "Dot Plot":
                fig = plot_gsea_dotplot(gsea_df, top_n_pathways)
            elif gsea_plot_type == "NES Bar Plot":
                fig = plot_gsea_barplot(gsea_df, top_n_pathways)
            else:
                fig = plot_gsea_waterfall(gsea_df, top_n_pathways)

            st.plotly_chart(fig, use_container_width=True)

            # GSEA results table
            with st.expander("📄 GSEA Results Table", expanded=False):
                st.dataframe(gsea_df, use_container_width=True, height=400)

            # Download
            csv_gsea = gsea_df.to_csv(index=False, sep='\t')
            st.download_button(
                "⬇️ Download GSEA Results (TSV)",
                csv_gsea, "kegg_gsea_results.tsv",
                "text/tab-separated-values", key="dl_gsea"
            )


# ═══════════════════════════════════════════════════════════════════════
#  TAB 4: COG ANALYSIS
# ═══════════════════════════════════════════════════════════════════════
with tab_cog:
    st.markdown('<div class="section-header">🏷️ COG Functional Category Analysis</div>',
                unsafe_allow_html=True)

    if not st.session_state.analysis_done:
        st.info("👈 Click **Run Full Analysis** in the sidebar to start.")
    else:
        cog_mapping = st.session_state.cog_mapping
        n_cog_mapped = len(cog_mapping)

        st.markdown(f"""
        <div class="info-panel">
            <b>COG Mapping:</b> {n_cog_mapped} genes with COG category assignments
            — Source: {cog_source}
        </div>
        """, unsafe_allow_html=True)

        if n_cog_mapped == 0:
            st.warning("No COG category assignments found. Try a different COG source.")
        else:
            # ── Distributions ──────────────────────────────────────────
            st.markdown("#### COG Category Distribution")

            all_gene_ids = df['gene_id'].dropna().unique().tolist()
            up_ids = up_genes['gene_id'].dropna().tolist()
            down_ids = down_genes['gene_id'].dropna().tolist()

            dist_all = get_cog_distribution(all_gene_ids, cog_mapping, "All Genes")
            dist_up = get_cog_distribution(up_ids, cog_mapping, "Upregulated")
            dist_down = get_cog_distribution(down_ids, cog_mapping, "Downregulated")

            cog_dist_plot = st.selectbox(
                "Distribution plot type",
                ["Grouped Bar Chart", "Stacked Bar Chart", "Pie Chart", "Heatmap"],
                key="cog_dist_type"
            )

            if cog_dist_plot == "Grouped Bar Chart":
                fig = plot_cog_distribution_bar([dist_all, dist_up, dist_down])
            elif cog_dist_plot == "Stacked Bar Chart":
                fig = plot_cog_distribution_stacked([dist_all, dist_up, dist_down])
            elif cog_dist_plot == "Pie Chart":
                pie_group = st.radio("Show pie for:", ["All Genes", "Upregulated", "Downregulated"],
                                      horizontal=True, key="cog_pie_group")
                if pie_group == "All Genes":
                    fig = plot_cog_pie(dist_all, f"COG Proportions — {pie_group}")
                elif pie_group == "Upregulated":
                    fig = plot_cog_pie(dist_up, f"COG Proportions — {pie_group}")
                else:
                    fig = plot_cog_pie(dist_down, f"COG Proportions — {pie_group}")
            else:
                fig = plot_cog_heatmap(dist_up, dist_down)

            st.plotly_chart(fig, use_container_width=True)

            # ── Enrichment ─────────────────────────────────────────────
            st.markdown("#### COG Category Enrichment")

            cog_enrich_dir = st.radio(
                "Enrichment for:",
                ["All DE genes", "Upregulated", "Downregulated"],
                horizontal=True, key="cog_enrich_dir"
            )

            if cog_enrich_dir == "All DE genes":
                cog_enrich_df = st.session_state.cog_enrich_all
            elif cog_enrich_dir == "Upregulated":
                cog_enrich_df = st.session_state.cog_enrich_up
            else:
                cog_enrich_df = st.session_state.cog_enrich_down

            if cog_enrich_df.empty:
                st.info("No enriched COG categories found for this selection.")
            else:
                cog_plot_type = st.selectbox(
                    "Enrichment plot type",
                    ["Dot Plot", "Bar Plot"],
                    key="cog_enrich_plot"
                )

                if cog_plot_type == "Dot Plot":
                    fig = plot_cog_enrichment_dotplot(cog_enrich_df, top_n_pathways,
                                                      f"COG Enrichment — {cog_enrich_dir}")
                else:
                    fig = plot_cog_enrichment_bar(cog_enrich_df, top_n_pathways,
                                                  f"COG Enrichment — {cog_enrich_dir}")

                st.plotly_chart(fig, use_container_width=True)

            # COG Category Reference Table
            with st.expander("📖 COG Category Reference", expanded=False):
                ref_data = []
                for cat, desc in COG_CATEGORIES.items():
                    # Find supercategory
                    supercat = "Other"
                    for sc, cats in COG_SUPERCATEGORIES.items():
                        if cat in cats:
                            supercat = sc
                            break
                    ref_data.append({
                        'Category': cat,
                        'Description': desc,
                        'Super-category': supercat,
                    })
                st.dataframe(pd.DataFrame(ref_data), use_container_width=True, height=400)

            # Enrichment results table
            with st.expander("📄 COG Enrichment Results Table", expanded=False):
                st.dataframe(cog_enrich_df, use_container_width=True, height=400)

            # Downloads
            col_dl1, col_dl2 = st.columns(2)
            with col_dl1:
                csv_dist = pd.concat([dist_all, dist_up, dist_down]).to_csv(index=False, sep='\t')
                st.download_button(
                    "⬇️ Download COG Distribution (TSV)",
                    csv_dist, "cog_distribution.tsv",
                    "text/tab-separated-values", key="dl_cog_dist"
                )
            with col_dl2:
                if not cog_enrich_df.empty:
                    csv_enrich = cog_enrich_df.to_csv(index=False, sep='\t')
                    st.download_button(
                        "⬇️ Download COG Enrichment (TSV)",
                        csv_enrich, f"cog_enrichment_{cog_enrich_dir.lower().replace(' ', '_')}.tsv",
                        "text/tab-separated-values", key="dl_cog_enrich"
                    )


# ═══════════════════════════════════════════════════════════════════════
#  TAB 5: COMBINED SUMMARY
# ═══════════════════════════════════════════════════════════════════════
with tab_combined:
    st.markdown('<div class="section-header">📊 Combined Analysis Summary</div>',
                unsafe_allow_html=True)

    if not st.session_state.analysis_done:
        st.info("👈 Click **Run Full Analysis** in the sidebar to start.")
    else:
        # Summary metrics
        st.markdown("#### Analysis Overview")

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            n_kegg = len(st.session_state.kegg_ora_all[
                st.session_state.kegg_ora_all['p.adjust'] < padj_cutoff
            ]) if not st.session_state.kegg_ora_all.empty else 0
            st.metric("Enriched KEGG Pathways", n_kegg)
        with col2:
            n_gsea = len(st.session_state.kegg_gsea_results) if not st.session_state.kegg_gsea_results.empty else 0
            st.metric("GSEA Pathways Tested", n_gsea)
        with col3:
            n_cog_sig = len(st.session_state.cog_enrich_all[
                st.session_state.cog_enrich_all['p.adjust'] < padj_cutoff
            ]) if not st.session_state.cog_enrich_all.empty else 0
            st.metric("Enriched COG Categories", n_cog_sig)
        with col4:
            st.metric("Mapped Genes", len(st.session_state.id_mapping))

        # Side-by-side comparison
        fig_combined = plot_kegg_cog_summary(
            st.session_state.kegg_ora_all,
            st.session_state.cog_enrich_all,
        )
        st.plotly_chart(fig_combined, use_container_width=True)

        # KEGG + GSEA combined view
        st.markdown("#### Top KEGG Results Comparison")

        col_k1, col_k2 = st.columns(2)
        with col_k1:
            st.markdown("**ORA — Upregulated**")
            if not st.session_state.kegg_ora_up.empty:
                display_df = st.session_state.kegg_ora_up[['Description', 'Count', 'p.adjust']].head(10)
                display_df['p.adjust'] = display_df['p.adjust'].apply(lambda x: f"{x:.2e}")
                st.dataframe(display_df, use_container_width=True, hide_index=True)
            else:
                st.info("No enriched pathways")

        with col_k2:
            st.markdown("**ORA — Downregulated**")
            if not st.session_state.kegg_ora_down.empty:
                display_df = st.session_state.kegg_ora_down[['Description', 'Count', 'p.adjust']].head(10)
                display_df['p.adjust'] = display_df['p.adjust'].apply(lambda x: f"{x:.2e}")
                st.dataframe(display_df, use_container_width=True, hide_index=True)
            else:
                st.info("No enriched pathways")

        # Download all results as Excel
        st.markdown("#### 📦 Download All Results")

        excel_buf = io.BytesIO()
        with pd.ExcelWriter(excel_buf, engine='openpyxl') as writer:
            # Processed data
            df.to_excel(writer, sheet_name='Processed_Data', index=False)
            # KEGG ORA
            if not st.session_state.kegg_ora_all.empty:
                st.session_state.kegg_ora_all.to_excel(writer, sheet_name='KEGG_ORA_All', index=False)
            if not st.session_state.kegg_ora_up.empty:
                st.session_state.kegg_ora_up.to_excel(writer, sheet_name='KEGG_ORA_Up', index=False)
            if not st.session_state.kegg_ora_down.empty:
                st.session_state.kegg_ora_down.to_excel(writer, sheet_name='KEGG_ORA_Down', index=False)
            # GSEA
            if not st.session_state.kegg_gsea_results.empty:
                st.session_state.kegg_gsea_results.to_excel(writer, sheet_name='KEGG_GSEA', index=False)
            # COG
            if not st.session_state.cog_enrich_all.empty:
                st.session_state.cog_enrich_all.to_excel(writer, sheet_name='COG_Enrichment_All', index=False)
            if not st.session_state.cog_enrich_up.empty:
                st.session_state.cog_enrich_up.to_excel(writer, sheet_name='COG_Enrichment_Up', index=False)
            if not st.session_state.cog_enrich_down.empty:
                st.session_state.cog_enrich_down.to_excel(writer, sheet_name='COG_Enrichment_Down', index=False)

        st.download_button(
            "📥 Download Complete Results (Excel)",
            excel_buf.getvalue(),
            "functional_enrichment_results.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True,
        )

        st.markdown("""
        <div class="info-panel">
            <b>Excel workbook contains:</b>
            <ul>
                <li><b>Processed_Data</b> — Parsed DESeq2 results with extracted gene info</li>
                <li><b>KEGG_ORA_*</b> — Over-Representation Analysis results (All/Up/Down)</li>
                <li><b>KEGG_GSEA</b> — Gene Set Enrichment Analysis results</li>
                <li><b>COG_Enrichment_*</b> — COG category enrichment results (All/Up/Down)</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════════
#  FOOTER
# ═══════════════════════════════════════════════════════════════════════
st.markdown("---")
st.markdown(
    "<div style='text-align:center; color:#aaa; font-size:0.8rem; padding:1rem;'>"
    "Functional Enrichment & GSEA App · Built with Streamlit · "
    "Analysis: KEGG REST API, Fisher's exact test, gseapy · "
    "Inspired by <a href='https://github.com/RickGelhausen/pathsnake'>pathsnake</a>"
    "</div>",
    unsafe_allow_html=True,
)
