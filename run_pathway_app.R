#Script to run the script once the app has been cloned from the github repository 

library(dplyr)
library(shiny)
library(tidyverse)
library(readxl)
library(ggplot2)

#replace your_path below with the path to the folder storing the cloned pathway_analysis github repository 
setwd("your_path/pathway_analysis")
getwd()

shiny::runApp("shiny_app.R")