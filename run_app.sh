#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Functional Enrichment & GSEA Analysis App — Launch Script
# ═══════════════════════════════════════════════════════════════

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ENV_NAME="enrichment_app"

# ── Check if conda env exists ─────────────────────────────────
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "╔══════════════════════════════════════════════════════╗"
    echo "║  Creating conda environment: ${ENV_NAME}            ║"
    echo "║  This may take a few minutes on first run...        ║"
    echo "╚══════════════════════════════════════════════════════╝"
    conda env create -f environment.yml
fi

# ── Activate and run ──────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║  Starting Functional Enrichment & GSEA App          ║"
echo "║  Open: http://localhost:8501                        ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""

# Use conda run to execute within the environment
conda run -n "$ENV_NAME" streamlit run app.py \
    --server.port 8501 \
    --server.headless true \
    --browser.gatherUsageStats false \
    --theme.primaryColor "#3498db" \
    --theme.backgroundColor "#ffffff" \
    --theme.secondaryBackgroundColor "#f0f2f6" \
    --theme.textColor "#1a1a2e"
