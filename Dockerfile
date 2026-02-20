# Base image with R and Shiny pre-installed (R 4.5+ required for limpa)
FROM rocker/shiny:4.5.0

# 1. Install System Dependencies (Linux libraries needed for R packages)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    cmake \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*

# 2. Install R Dependencies (CRAN)
RUN R -e "install.packages(c('bslib', 'readr', 'tibble', 'dplyr', 'tidyr', 'ggplot2', 'httr2', 'rhandsontable', 'DT', 'arrow', 'shinyjs', 'plotly', 'stringr', 'ggrepel', 'remotes', 'BiocManager', 'markdown', 'shinyFiles', 'jsonlite'), repos='https://cloud.r-project.org/')"

# 2b. Install phosphoproteomics packages (KSEA kinase activity, sequence logos)
RUN R -e "install.packages(c('KSEAapp', 'ggseqlogo'), repos='https://cloud.r-project.org/')"

# 2c. Install graphics/font dependencies for clusterProfiler/enrichplot
RUN R -e "install.packages(c('systemfonts', 'gdtools', 'Rcpp'), repos='https://cloud.r-project.org/')"

# 2d. Install network visualization dependencies for enrichplot
RUN R -e "install.packages(c('ggraph', 'graphlayouts', 'tidygraph', 'scatterpie', 'shadowtext', 'ggforce'), repos='https://cloud.r-project.org/')"

# 3. Install Bioconductor Packages (Specific versions for stability)
# Install in correct dependency order to avoid compilation issues
RUN R -e "BiocManager::install(c('DOSE', 'GOSemSim', 'yulab.utils'), ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install(c('limma', 'limpa', 'ComplexHeatmap', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db', 'ggridges'), ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install(c('ggtree', 'ggtangle'), ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install(c('clusterProfiler', 'enrichplot'), ask=FALSE, update=FALSE)"

# 4. Copy the App Files into the image
COPY app.R /srv/shiny-server/app.R
COPY R/ /srv/shiny-server/R/
COPY contaminants/ /srv/shiny-server/contaminants/

# (Optional: If you have your logo file locally, uncomment the next line)
# COPY funny_scientist.png /srv/shiny-server/funny_scientist.png

# 5. Expose the port (Shiny runs on 3838 by default)
EXPOSE 3838

# 6. Run the App
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 7860)"]
