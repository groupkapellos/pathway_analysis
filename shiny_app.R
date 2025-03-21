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
             ),
             fluidRow(
               column(12,
                      h3("Gene and Pathway Database"),
                      tableOutput("gene_pathway_table")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  # Handle the pathway search (Search by Pathway)
  observeEvent(input$search, {
    req(input$pathway)
    
    if (input$pathway == "Select a Pathway") {
      output$results <- renderTable({
        data.frame(Message = "Please select a pathway to search.")
      })
    } else {
      filtered_pathway_data <- data %>%
        filter(pathway == input$pathway) %>%
        select(pathway, Positive_Genes)
      
      output$results <- renderTable({
        filtered_pathway_data
      })
    }
  })
  
  filtered_data <- reactiveVal()
  
  observeEvent(input$find_pathways, {
    req(input$gene_list)
    
    gene_list <- unlist(strsplit(input$gene_list, ","))
    gene_list <- trimws(gene_list)
    
    unmatched_genes <- gene_list
    genes_not_associated <- c()  
    
    for (gene in gene_list) {
      gene_found <- any(sapply(data$Positive_Genes, function(positive_genes) {
        gene %in% trimws(strsplit(positive_genes, ",")[[1]])
      }))
      
      if (!gene_found) {
        genes_not_associated <- c(genes_not_associated, gene)
      }
    }
    
    unmatched_genes_sorted <- sort(genes_not_associated)
    
    output$unmatched_message <- renderText({
      if (length(unmatched_genes_sorted) > 0) {
        paste("The following genes were not found to be associated with the database: ", 
              paste(unmatched_genes_sorted, collapse = ", "))
      } else {
        "All entered genes were found to be associated with pathways in the database."
      }
    })
    
    matched_pathways <- data %>%
      filter(sapply(strsplit(Positive_Genes, ", "), function(genes) {
        any(gene_list %in% genes)
      })) %>%
      select(pathway, Positive_Genes)
    
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Positive_Genes, function(genes) {
      intersect(gene_list, strsplit(genes, ", ")[[1]])
    })
    
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Detected_Genes, function(genes) {
      if (length(genes) > 0) {
        return(paste(sort(genes), collapse = ", "))
      } else {
        return(NULL)
      }
    })
    
    matched_pathways$Detected_Gene_Count <- sapply(matched_pathways$Detected_Genes, function(detected_genes) {
      if (!is.null(detected_genes) && is.character(detected_genes)) {
        return(length(strsplit(detected_genes, ", ")[[1]]))
      } else {
        return(0)
      }
    })
    
    # Frequency Score Calculation
    pathway_counts <- data %>%
      group_by(pathway) %>%
      summarise(
        total_positive = sum(!is.na(Positive_Genes)),  # Count non-NA values as positives
        total_negative = sum(is.na(Positive_Genes)),   # Count NA values as negatives
        .groups = "drop"
      ) %>%
      mutate(Frequency_Score = total_positive / (total_positive + total_negative))
    
    matched_pathways <- matched_pathways %>%
      left_join(pathway_counts, by = "pathway") %>%
      arrange(desc(Detected_Gene_Count))
    
    filtered_data(matched_pathways)
    
    pathway_gene_count <- matched_pathways %>%
      mutate(gene_count = sapply(strsplit(Detected_Genes, ", "), length)) %>%
      group_by(pathway) %>%
      summarise(total_genes = sum(gene_count), Frequency_Score = mean(Frequency_Score, na.rm = TRUE)) %>%
      arrange(desc(total_genes)) %>%
      top_n(10, total_genes)
    
    output$gene_count_plot <- renderPlot({
      ggplot(pathway_gene_count, aes(x = reorder(pathway, -total_genes), y = total_genes)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        geom_text(aes(label = total_genes), vjust = -0.3, size = 5, color = "black") +
        labs(x = "Pathway", y = "Number of Genes", title = "Top 10 Pathways with Associated Genes") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5)
        ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
    }, height = 400)
    
    output$pathway_results <- renderTable({
      matched_pathways %>%
        filter(!is.null(Detected_Genes)) %>%
        select(pathway, Detected_Genes, Detected_Gene_Count, Frequency_Score)
    })
    
    upset_data <- matched_pathways %>%
      mutate(pathway_gene_set = strsplit(Detected_Genes, ", ")) %>%
      select(pathway, pathway_gene_set) %>%
      unnest(pathway_gene_set) %>%
      group_by(pathway, pathway_gene_set) %>%
      summarise(count = n(), .groups = 'drop')
    
    upset_matrix <- upset_data %>%
      pivot_wider(names_from = pathway, values_from = count, values_fill = list(count = 0)) %>%
      select(-pathway_gene_set)
    
    output$upset_plot <- renderPlot({
      upset(
        as.data.frame(upset_matrix), 
        sets.bar.color = "skyblue", 
        main.bar.color = "darkblue",
        order.by = "freq", 
        set_size.show = TRUE,
        text.scale = 1.5
      )
    })
    
    gene_pathway_data <- matched_pathways %>%
      separate_rows(Detected_Genes, sep = ", ") %>%
      select(Detected_Genes, pathway) %>%
      distinct() %>%
      group_by(Detected_Genes) %>%
      summarise(pathways = paste(pathway, collapse = ", "), .groups = "drop")
    
    output$gene_pathway_table <- renderTable({
      gene_pathway_data
    })
  })
  
  output$download_filtered_data <- downloadHandler(
    filename = function() {
      paste("filtered_pathway_gene_data.csv")
    },
    content = function(file) {
      write.csv(filtered_data(), file)
    }
  )
}




# Run the Shiny application
shinyApp(ui = ui, server = server)
