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
                      # Dropdown for pathway selection, sorted alphabetically
                      selectInput("pathway", "Select Pathway Name:", 
                                  choices = c("Select a Pathway", sort(unique(data$pathway))), 
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
                      h3("Gene Search Results"),
                      plotOutput("gene_count_plot"),    # Barplot output - now first
                      tableOutput("pathway_results"),   # Table output for pathways
                      textOutput("unmatched_message")  # Text message for unmatched genes
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  # Handle the gene list search
  observeEvent(input$find_pathways, {
    req(input$gene_list)
    
    # Split the input gene list into individual genes
    gene_list <- unlist(strsplit(input$gene_list, ","))
    gene_list <- trimws(gene_list)  # Remove any leading/trailing spaces
    
    # Track genes not found in any pathway
    unmatched_genes <- gene_list
    genes_not_associated <- c()  # To track genes that are not associated with any pathway
    
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
    
    # Sort the genes in Detected_Genes column alphabetically for each row
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Detected_Genes, function(genes) {
      if (length(genes) > 0) {
        return(paste(sort(genes), collapse = ", "))  # Sort genes alphabetically for each pathway
      } else {
        return(NULL)  # Return NULL instead of "No detected genes"
      }
    })
    
    # Remove detected genes from unmatched_genes list
    detected_genes <- unique(unlist(matched_pathways$Detected_Genes))
    unmatched_genes <- setdiff(gene_list, detected_genes)
    
    # Track genes that were not associated with any pathway
    for (gene in gene_list) {
      if (!any(gene %in% detected_genes)) {
        genes_not_associated <- c(genes_not_associated, gene)
      }
    }
    
    # Fix the mutate issue: check if Detected_Genes is character before using strsplit
    matched_pathways$Detected_Gene_Count <- sapply(matched_pathways$Detected_Genes, function(detected_genes) {
      if (!is.null(detected_genes) && is.character(detected_genes)) {
        return(length(strsplit(detected_genes, ", ")[[1]]))  # Ensure count is a whole number
      } else {
        return(0)
      }
    })
    
    # Sort by Detected_Gene_Count in descending order
    matched_pathways <- matched_pathways %>%
      arrange(desc(Detected_Gene_Count))
    
    # Limit the plot to the top 10 pathways based on Detected_Gene_Count, excluding "Not Found"
    pathway_gene_count <- matched_pathways %>%
      mutate(gene_count = sapply(strsplit(Detected_Genes, ", "), length)) %>%
      group_by(pathway) %>%
      summarise(total_genes = sum(gene_count)) %>%
      arrange(desc(total_genes)) %>%
      top_n(10, total_genes)  # Limit to the top 10 pathways
    
    # Output the bar plot with top 10 pathways ordered by total gene count
    output$gene_count_plot <- renderPlot({
      ggplot(pathway_gene_count, aes(x = reorder(pathway, -total_genes), y = total_genes)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        labs(x = "Pathway", y = "Number of Genes", title = "Top 10 Pathways with Associated Genes") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase size of x-axis text
          axis.text.y = element_text(size = 14),  # Increase size of y-axis text
          axis.title.x = element_text(size = 16),  # Increase size of x-axis title
          axis.title.y = element_text(size = 16),  # Increase size of y-axis title
          plot.title = element_text(size = 18, hjust = 0.5)  # Increase size of plot title
        )
    })
    
    # Display the unmatched genes as a message
    output$unmatched_message <- renderText({
      if (length(unmatched_genes) > 0) {
        # Sort the unmatched genes alphabetically
        unmatched_genes_sorted <- sort(unmatched_genes)
        paste("The following genes were not found to be associated with the database: ", 
              paste(unmatched_genes_sorted, collapse = ", "))
      } else {
        "All entered genes were found to be associated with pathways in the database."
      }
    })
    
    # Output the results of pathways associated with the gene list (excluding "Not Found" row)
    output$pathway_results <- renderTable({
      matched_pathways %>%
        filter(!is.null(Detected_Genes)) %>%  # Remove rows where Detected_Genes is NULL
        select(pathway, Detected_Genes, Detected_Gene_Count)  # Display relevant columns
    })
    
  })
}

# Run the application
shinyApp(ui = ui, server = server)


