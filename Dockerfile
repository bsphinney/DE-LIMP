# DE-LIMP Application Image
# Uses pre-built base image with all R/Bioconductor/MOFA2 dependencies
# Base image: Dockerfile.base (rebuild only when dependencies change)
FROM brettphinney/delimp-base:v3.1

# Install clusterProfiler + enrichplot + org.Hs.eg.db (required for GSEA)
# Cache-bust: v2
RUN R -e "options(repos = BiocManager::repositories()); install.packages(c('clusterProfiler','enrichplot','org.Hs.eg.db'), quiet=TRUE)" 2>&1 | tail -5 || true

# Install ggdendro for Data Completeness dendrogram visualization
RUN R -e "if (!requireNamespace('ggdendro', quietly=TRUE)) install.packages('ggdendro', repos='https://cloud.r-project.org/')" 2>/dev/null || true

# Copy the App Files into the image
COPY app.R /srv/shiny-server/app.R
COPY R/ /srv/shiny-server/R/
COPY contaminants/ /srv/shiny-server/contaminants/
COPY VERSION /srv/shiny-server/VERSION
COPY stats/ /srv/shiny-server/stats/

# Expose the port
EXPOSE 3838

# Run the App
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 7860)"]
