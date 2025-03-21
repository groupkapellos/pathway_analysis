# Load necessary libraries
library(shiny)
library(dplyr)
library(readr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(UpSetR)
library(stringr)

# Read the Excel file
data <- read_excel("data/compiled_gene_lists.xlsx")

ui <- fluidPage(
  titlePanel("Pathway-Gene Association Mapper"),
  
  # Definition Box
  wellPanel(
    h3("Gene Definitions"),
    p(strong("Positive Genes:"), " Database-derived genes that are truly associated with a pathway or biological process."),
    p(strong("Negative Genes:"), " Database-derived genes that are only loosely associated with a pathway or biological process.")
  ),
  
  tabsetPanel(
    tabPanel("Search by Pathway",
             fluidRow(
               column(12,
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
    tabPanel("Search Gene List",
             fluidRow(
               column(12,
                      textAreaInput("gene_list", "Enter Genes (comma-separated):", ""),
                      actionButton("find_pathways", "Find Pathways")
               ),
               column(12,
                      h3("Gene Search Results"),
                      plotOutput("gene_count_plot"), 
                      tableOutput("pathway_results"),  
                      textOutput("unmatched_message") 
               ),
               column(12,
                      downloadButton("download_filtered_data", "Download Filtered Data")
               )
             ),
             fluidRow(
               column(12,
                      h3("UpSet Plot - Gene and Pathway Association"),
                      plotOutput("upset_plot")
               )
             )
    )
  )
)


server <- function(input, output, session) {
  
  # Handle the pathway search (Search by Pathway)
  observeEvent(input$search, {
    req(input$pathway)
    
    # If no pathway is selected, display an empty table
    if (input$pathway == "Select a Pathway") {
      output$results <- renderTable({
        data.frame(Message = "Please select a pathway to search.")
      })
    } else {
      # Filter the data for the selected pathway
      filtered_pathway_data <- data %>%
        filter(pathway == input$pathway) %>%
        select(pathway, Positive_Genes)
      
      #cam we change the output names to pathway and positive genes 
      # Show the filtered results in the results table
      output$results <- renderTable({
        filtered_pathway_data
      })
    }
  })
  
  # Store the filtered data in a reactive value
  filtered_data <- reactiveVal()
  
  # Handle the gene list search
  observeEvent(input$find_pathways, {
    req(input$gene_list)
    
    # Split the input gene list into individual genes and clean up white spaces
    gene_list <- unlist(strsplit(input$gene_list, ","))
    gene_list <- trimws(gene_list)  # Trim spaces
    
    # Track genes not found in any pathway
    unmatched_genes <- gene_list
    genes_not_associated <- c()  
    
    # Check each gene in the input list to see if it is found in the Positive_Genes column
    for (gene in gene_list) {
      # Check if the gene is present in any of the pathways' Positive_Genes
      gene_found <- any(sapply(data$Positive_Genes, function(positive_genes) {
        gene %in% trimws(strsplit(positive_genes, ",")[[1]])  # Check if gene is in Positive_Genes
      }))
      
      # If the gene was not found, add it to the unmatched list
      if (!gene_found) {
        genes_not_associated <- c(genes_not_associated, gene)
      }
    }
    
    # Sort unmatched genes alphabetically
    unmatched_genes_sorted <- sort(genes_not_associated)
    
    # Display the unmatched genes as a message
    output$unmatched_message <- renderText({
      if (length(unmatched_genes_sorted) > 0) {
        paste("The following genes were not found to be associated with the database: ", 
              paste(unmatched_genes_sorted, collapse = ", "))
      } else {
        "All entered genes were found to be associated with pathways in the database."
      }
    })
    
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
    
    # Store the matched pathways data into a reactive value (so it can be downloaded)
    filtered_data(matched_pathways)
    
    # Limit the plot to the top 10 pathways based on Detected_Gene_Count
    pathway_gene_count <- matched_pathways %>%
      mutate(gene_count = sapply(strsplit(Detected_Genes, ", "), length)) %>%
      group_by(pathway) %>%
      summarise(total_genes = sum(gene_count)) %>%
      arrange(desc(total_genes)) %>%
      top_n(10, total_genes)  # Limit to the top 10 pathways
    
    output$gene_count_plot <- renderPlot({
      ggplot(pathway_gene_count, aes(x = reorder(pathway, -total_genes), y = total_genes)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        geom_text(aes(label = total_genes), vjust = -0.3, size = 4, color = "black") +  # Add numbers to bars
        labs(x = "Pathway", y = "Number of Genes", title = "Top 10 Pathways with Associated Genes") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5),
          plot.margin = margin(t = 20, r = 20, b = 30, l = 20),
          panel.background = element_blank(),  # Remove the background
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()   # Remove minor grid lines
        )
    })
    
    # Output the results of pathways associated with the gene list (excluding "Not Found" row)
    output$pathway_results <- renderTable({
      matched_pathways %>%
        filter(!is.null(Detected_Genes)) %>%
        select(pathway, Detected_Genes, Detected_Gene_Count)
    })
    
    # Prepare the data for the UpSet plot
    upset_data <- matched_pathways %>%
      mutate(pathway_gene_set = strsplit(Detected_Genes, ", ")) %>%
      select(pathway, pathway_gene_set) %>%
      unnest(pathway_gene_set) %>%
      group_by(pathway, pathway_gene_set) %>%
      summarise(count = n(), .groups = 'drop')  # Ensure summarise works correctly
    
    # Pivot the data to create the binary matrix for the UpSet plot using pivot_wider
    upset_matrix <- upset_data %>%
      pivot_wider(names_from = pathway, values_from = count, values_fill = list(count = 0)) %>%
      select(-pathway_gene_set)  # Remove unnecessary columns
    
    # Debugging step: Print the transformed upset matrix to the console
    print(upset_matrix)
    
    # Render the UpSet plot with all pathways included
    output$upset_plot <- renderPlot({
      upset(as.data.frame(upset_matrix), 
            sets.bar.color = "skyblue", 
            main.bar.color = "darkblue",
            order.by = "freq", 
            set_size.show = TRUE)
      
    })
  })
  
  # Set up the download handler for the filtered data
  output$download_filtered_data <- downloadHandler(
    filename = function() {
      paste("filtered_pathway_gene_data.csv")
    },
    content = function(file) {
      write.csv(filtered_data(), file)
    }
  )
}


#but there is still a problem in the uspet plot. Always the second row from the filtered data is not appearing 


shinyApp(ui = ui, server = server)

