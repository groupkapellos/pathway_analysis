# Load necessary libraries
library(shiny)
library(dplyr)
library(readr)
library(tidyverse)
library(readxl)
library(ggplot2)

# Read the Excel file
data <- read_excel("./data/compiled_gene_lists.xlsx")

# Define UI for the application
ui <- fluidPage(
  titlePanel("Pathway Gene Finder"),
  
  # Tabs for both functionalities (Pathway search and Gene list search)
  tabsetPanel(
    tabPanel("Pathway Search",
             fluidRow(
               column(12,
                      # Dropdown for pathway selection
                      selectInput("pathway", "Select Pathway Name:", 
                                  choices = c("Select a Pathway", unique(data$pathway)), 
                                  selected = "Select a Pathway"),
                      actionButton("search", "Search")
               ),
               column(12,
                      h3("Search Results"),
                      tableOutput("results")
               )
             )
    ),
    tabPanel("Gene List Search",
             fluidRow(
               column(12,
                      # Text input or file input for user to enter a gene list
                      textAreaInput("gene_list", "Enter Genes (comma-separated):", ""),
                      actionButton("find_pathways", "Find Pathways")
               ),
               column(12,
                      h3("Associated Pathways"),
                      tableOutput("pathway_results"),
                      plotOutput("gene_count_plot") # Barplot output
               )
             )
    )
  )
)

# Define server logic for the application
server <- function(input, output, session) {
  
  # Handle the gene list search
  observeEvent(input$find_pathways, {
    req(input$gene_list)
    
    # Split the input gene list into individual genes
    gene_list <- unlist(strsplit(input$gene_list, ","))
    gene_list <- trimws(gene_list)  # Remove any leading/trailing spaces
    
    # Filter the data for pathways that contain any of the genes in the gene list
    matched_pathways <- data %>%
      filter(sapply(strsplit(Positive_Genes, ", "), function(genes) {
        any(gene_list %in% genes)
      })) %>%
      select(pathway, Positive_Genes)  # Select only relevant columns
    
    # Show only the genes that were actually entered (i.e., that exist in the Positive_Genes slot)
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Positive_Genes, function(genes) {
      intersect(gene_list, strsplit(genes, ", ")[[1]])
    })
    
    # Ensure that Detected_Genes is converted into a proper format (if it's a list, turn it into a string)
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Detected_Genes, function(genes) {
      if (length(genes) > 0) {
        return(paste(genes, collapse = ", "))
      } else {
        return("No detected genes")
      }
    })
    
    # Count how many genes were detected for each pathway (convert the count to numeric)
    matched_pathways$Detected_Gene_Count <- sapply(matched_pathways$Detected_Genes, function(detected_genes) {
      length(strsplit(detected_genes, ", ")[[1]])
    })
    
    # Add a frequency score: Number of detected genes vs total genes for that pathway
    matched_pathways$Frequency_Score <- sapply(matched_pathways$Positive_Genes, function(genes) {
      detected_genes <- intersect(gene_list, strsplit(genes, ", ")[[1]])
      detected_count <- length(detected_genes)
      total_count <- length(strsplit(genes, ", ")[[1]])
      detected_count / total_count  # Frequency score as ratio
    })
    
    # Order pathways by the number of detected genes in descending order
    matched_pathways <- matched_pathways %>%
      arrange(desc(Detected_Gene_Count))
    
    # Debugging: print structure of matched_pathways to check if it’s a data frame
    print(str(matched_pathways))
    
    # Ensure that it's a data frame before rendering
    matched_pathways <- as.data.frame(matched_pathways)
    
    # Output the results of pathways associated with the gene list
    output$pathway_results <- renderTable({
      matched_pathways %>%
        select(pathway, Detected_Genes, Detected_Gene_Count, Frequency_Score)  # Display relevant columns
    })
    
    # Create a barplot showing the number of genes in each pathway
    pathway_gene_count <- matched_pathways %>%
      mutate(gene_count = sapply(strsplit(Positive_Genes, ", "), length)) %>%
      group_by(pathway) %>%
      summarise(total_genes = sum(gene_count)) %>%
      arrange(desc(total_genes))
    
    output$gene_count_plot <- renderPlot({
      ggplot(pathway_gene_count, aes(x = reorder(pathway, -total_genes), y = total_genes)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        labs(x = "Pathway", y = "Number of Genes", title = "Pathways with Associated Genes") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase size of x-axis text
          axis.text.y = element_text(size = 14),  # Increase size of y-axis text
          axis.title.x = element_text(size = 16),  # Increase size of x-axis title
          axis.title.y = element_text(size = 16),  # Increase size of y-axis title
          plot.title = element_text(size = 18, hjust = 0.5)  # Increase size of plot title
        )
    })
    
  })
}

# Run the application
shinyApp(ui = ui, server = server)


#This also needs to work with multiple inputs, so the inputs are a comma separated gene list 
