# Use the official R base image from Docker Hub
FROM r-base

# Copy the R Shiny app files to the container
COPY src/app.R/ app/
COPY src/modules/ app/modules/
COPY requirements.txt app/

# Install necessary dependencies
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev

# Install required R packages
RUN R -e "install.packages(readLines('app/requirements.txt'), repos='https://cran.rstudio.com/')"

# Expose the Shiny app port
EXPOSE 3838

# Run the Shiny app on container startup
CMD ["R", "-e", "shiny::runApp('app/app.R', host = '0.0.0.0', port = 3838)"]