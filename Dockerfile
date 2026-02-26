# DE-LIMP Application Image
# Uses pre-built base image with all R/Bioconductor/MOFA2 dependencies
# Base image: Dockerfile.base (rebuild only when dependencies change)
FROM brettphinney/delimp-base:v3.1

# Copy the App Files into the image
COPY app.R /srv/shiny-server/app.R
COPY R/ /srv/shiny-server/R/
COPY contaminants/ /srv/shiny-server/contaminants/

# Expose the port
EXPOSE 3838

# Run the App
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 7860)"]
