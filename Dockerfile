# DE-LIMP Application Image
# Uses pre-built base image with all R/Bioconductor/MOFA2 dependencies
# Base image: Dockerfile.base (rebuild only when dependencies change)
FROM brettphinney/delimp-base:v3.1

# Install system deps for enrichplot (needs libmagick for magick R package)
# Cache-bust: v3
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    libmagick++-dev libharfbuzz-dev libfribidi-dev 2>/dev/null && \
    rm -rf /var/lib/apt/lists/*

# Install clusterProfiler + enrichplot (required for GSEA)
# Cache-bust: v6
RUN R -e " \
  cat('=== Installing ggtangle from GitHub ===\n'); \
  if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org/'); \
  remotes::install_github('YuLab-SMU/ggtangle', upgrade='never', force=TRUE); \
  cat('ggtangle installed:', requireNamespace('ggtangle', quietly=TRUE), '\n'); \
  cat('=== Installing enrichplot + clusterProfiler ===\n'); \
  options(repos = BiocManager::repositories()); \
  BiocManager::install(c('enrichplot','clusterProfiler'), ask=FALSE, update=FALSE); \
  cat('enrichplot installed:', requireNamespace('enrichplot', quietly=TRUE), '\n'); \
  cat('clusterProfiler installed:', requireNamespace('clusterProfiler', quietly=TRUE), '\n')" 2>&1 || true

# Install ggdendro for Data Completeness dendrogram visualization
RUN R -e "if (!requireNamespace('ggdendro', quietly=TRUE)) install.packages('ggdendro', repos='https://cloud.r-project.org/')" 2>/dev/null || true

# Copy the App Files into the image
COPY app.R /srv/shiny-server/app.R
COPY R/ /srv/shiny-server/R/
COPY contaminants/ /srv/shiny-server/contaminants/
COPY VERSION /srv/shiny-server/VERSION
COPY stats/ /srv/shiny-server/stats/
COPY config.yml /srv/shiny-server/config.yml

# Expose the port
EXPOSE 3838

# Run the App
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 7860)"]
