# Use the official R base image from Docker Hub
FROM r-base

# Copy the R script into the container
COPY ./ /app/

# Set the working directory
WORKDIR /app/src/multi_omics_platform

# Install required R packages
RUN Rscript -e "install.packages(c('magrittr','data.table'))"

# Set the ENTRYPOINT to handle the script execution
ENTRYPOINT ["Rscript", "main.R"]

# Set the CMD to provide the default argument value
CMD [filename]