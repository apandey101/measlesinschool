# Start from a basic devcontainers Ubuntu image
FROM mcr.microsoft.com/devcontainers/base:ubuntu

# Install R
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      r-base \
      r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# Install R packages you'll use
RUN R -e 'install.packages(c("dplyr", "Rcpp", "tidyr","ggplot2","plotly","shiny","tidyr","shinycssloaders","readxl","shinyBS","languageserver"), repos = "https://cloud.r-project.org")'

# Set default user (this matches what devcontainer.json expects)
USER vscode
