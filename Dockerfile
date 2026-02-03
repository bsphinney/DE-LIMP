# Base image with R and Shiny pre-installed
FROM rocker/shiny:4.2.1

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
    cmake \
    && rm -rf /var/lib/apt/lists/*

# 2. Install R Dependencies (CRAN)
RUN R -e "install.packages(c('bslib', 'readr', 'tibble', 'dplyr', 'tidyr', 'ggplot2', 'httr2', 'rhandsontable', 'DT', 'arrow', 'shinyjs', 'plotly', 'stringr', 'ggrepel', 'remotes', 'BiocManager'), repos='https://cloud.r-project.org/')"

# 3. Install Bioconductor Packages (Specific versions for stability)
RUN R -e "BiocManager::install(c('limma', 'limpa', 'ComplexHeatmap', 'clusterProfiler', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db', 'enrichplot', 'ggridges'), ask=FALSE)"

# 4. Copy the App Files into the image
COPY LIMP-D.R /srv/shiny-server/app.R
# (Optional: Copy your logo if you have one)
# COPY funny_scientist.png /srv/shiny-server/funny_scientist.png

# 5. Expose the port (Shiny runs on 3838 by default)
EXPOSE 3838

# 6. Run the App
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app.R', host = '0.0.0.0', port = 7860)"]
