# install_dependencies.R for installing packages needed for the pathway project app 

# List of required packages for the app
required_packages <- c(
  "shiny",      # For building the Shiny app
  "dplyr",      # For data manipulation
  "tidyverse",  # Includes several useful packages like ggplot2, tidyr, and readr
  "readxl",     # For reading Excel files
  "ggplot2"     # For creating visualizations (bar plot)
)

# Check which packages are already installed
installed_packages <- installed.packages()[,"Package"]

# Identify packages that are not installed
missing_packages <- required_packages[!(required_packages %in% installed_packages)]

# Install any missing packages
if(length(missing_packages)) {
  install.packages(missing_packages)
} else {
  message("All required packages are already installed.")
}

# Optional: Load all the required libraries after installation
invisible(lapply(required_packages, library, character.only = TRUE))

message("Dependencies installation complete!")
