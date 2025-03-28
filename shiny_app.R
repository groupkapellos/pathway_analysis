
# Load necessary libraries
library(shiny)
library(dplyr)
library(readr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(UpSetR)
library(stringr)

# Read the Excel file and filter for show == TRUE
data <- read_excel("data/pathway_gene_file.xlsx") %>%
  filter(!is.na(`Show`) & as.logical(`Show`) == TRUE)

ui <- fluidPage(
  titlePanel("Pathway-Gene Association Mapper"),
  
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
                      tableOutput("results"),
                      downloadButton("download_pathway_results", "Download Pathway Search Results")
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
                      tableOutput("gene_pathway_table"),
                      downloadButton("download_gene_pathway_table", "Download Gene-Pathway Table")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
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
    
    gene_list <- unlist(strsplit(input$gene_list, ",")) %>%
      trimws() %>%
      tolower()
    
    all_positive_pairs <- data %>%
      filter(!is.na(Positive_Genes)) %>%
      separate_rows(Positive_Genes, sep = ", ") %>%
      rename(Gene = Positive_Genes) %>%
      mutate(Gene = tolower(Gene))
    
    gene_pathway_matches <- all_positive_pairs %>%
      filter(Gene %in% gene_list) %>%
      group_by(Gene) %>%
      summarise(Matched_Pathways = n_distinct(pathway), .groups = "drop")
    
    gene_total_pathways <- all_positive_pairs %>%
      group_by(Gene) %>%
      summarise(Total_Pathways = n_distinct(pathway), .groups = "drop")
    
    gene_scores <- left_join(gene_pathway_matches, gene_total_pathways, by = "Gene") %>%
      mutate(
        Gene_Frequency_Score = Matched_Pathways / Total_Pathways,
        Match_Overview = paste0(Matched_Pathways, "/", Total_Pathways),
        Gene = toupper(Gene)
      )
    
    unmatched_genes <- setdiff(gene_list, unique(all_positive_pairs$Gene))
    
    output$unmatched_message <- renderText({
      if (length(unmatched_genes) > 0) {
        paste("The following genes were not found to be associated with the database:", 
              paste(toupper(unmatched_genes), collapse = ", "))
      } else {
        "All entered genes were found to be associated with pathways in the database."
      }
    })
    
    matched_pathways <- data %>%
      filter(sapply(strsplit(Positive_Genes, ", "), function(genes) {
        any(tolower(gene_list) %in% tolower(genes))
      })) %>%
      select(pathway, Positive_Genes)
    
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Positive_Genes, function(genes) {
      intersect(gene_list, tolower(strsplit(genes, ", ")[[1]]))
    })
    
    matched_pathways$Detected_Genes <- sapply(matched_pathways$Detected_Genes, function(genes) {
      if (length(genes) > 0) {
        return(toupper(paste(sort(genes), collapse = ", ")))
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
    
    pathway_gene_totals <- all_positive_pairs %>%
      group_by(pathway) %>%
      summarise(Total_Genes = n_distinct(Gene), .groups = "drop")
    
    matched_pathways$Detected_Gene_List <- sapply(matched_pathways$Detected_Genes, function(x) {
      strsplit(x, ", ")[[1]]
    })
    
    matched_pathways$Matched_Count <- sapply(matched_pathways$Detected_Gene_List, length)
    
    matched_pathways <- matched_pathways %>%
      left_join(pathway_gene_totals, by = "pathway") %>%
      mutate(
        Match_Overview = paste0(Matched_Count, "/", Total_Genes),
        Pathway_Frequency_Score = Matched_Count / Total_Genes
      )
    
    filtered_data(matched_pathways)
    
    pathway_gene_count <- matched_pathways %>%
      mutate(gene_count = sapply(strsplit(Detected_Genes, ", "), length)) %>%
      group_by(pathway) %>%
      summarise(total_genes = sum(gene_count), 
                Pathway_Frequency_Score = mean(Pathway_Frequency_Score, na.rm = TRUE)) %>%
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
        arrange(desc(Detected_Gene_Count)) %>%
        select(pathway, Detected_Genes, Detected_Gene_Count, Match_Overview, Pathway_Frequency_Score)
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
      top_pathways <- colSums(as.data.frame(upset_matrix))
      top_10_pathways <- names(sort(top_pathways, decreasing = TRUE))[1:min(10, length(top_pathways))]
      
      upset(
        as.data.frame(upset_matrix),
        sets = top_10_pathways,
        sets.bar.color = "skyblue", 
        main.bar.color = "darkblue",
        order.by = "freq", 
        set_size.show = TRUE,
        text.scale = 1.5,
        mainbar.y.label = "Top 10 Pathway Overlaps"
      )
    })
    
    gene_pathway_data <- matched_pathways %>%
      separate_rows(Detected_Genes, sep = ", ") %>%
      mutate(Detected_Genes = toupper(Detected_Genes)) %>%
      select(Detected_Genes, pathway) %>%
      distinct() %>%
      group_by(Detected_Genes) %>%
      summarise(pathways = paste(pathway, collapse = ", "), .groups = "drop")
    
    output$gene_pathway_table <- renderTable({
      gene_pathway_data %>%
        mutate(Pathway_Count = as.integer(str_count(pathways, ",") + 1)) %>%
        arrange(desc(Pathway_Count)) %>%
        select(Detected_Genes, pathways, Pathway_Count)
    })
    
    output$gene_frequency_table <- renderTable({
      gene_scores %>%
        select(Gene, Match_Overview, Gene_Frequency_Score)
    })
  })
  
  output$download_filtered_data <- downloadHandler(
    filename = function() {
      paste("filtered_pathway_gene_data.csv")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
  
  output$download_pathway_results <- downloadHandler(
    filename = function() {
      paste("pathway_search_results.csv")
    },
    content = function(file) {
      if (input$pathway != "Select a Pathway") {
        write.csv(data %>%
                    filter(pathway == input$pathway) %>%
                    select(pathway, Positive_Genes), file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No pathway selected."), file, row.names = FALSE)
      }
    }
  )
  
  output$download_gene_pathway_table <- downloadHandler(
    filename = function() {
      paste("gene_pathway_database.csv")
    },
    content = function(file) {
      gene_pathway_data <- filtered_data() %>%
        separate_rows(Detected_Genes, sep = ", ") %>%
        mutate(Detected_Genes = toupper(Detected_Genes)) %>%
        select(Detected_Genes, pathway) %>%
        distinct() %>%
        group_by(Detected_Genes) %>%
        summarise(pathways = paste(pathway, collapse = ", "), .groups = "drop") %>%
        mutate(Pathway_Count = as.integer(str_count(pathways, ",") + 1)) %>%
        arrange(desc(Pathway_Count)) %>%
        select(Detected_Genes, pathways, Pathway_Count)
      
      write.csv(gene_pathway_data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
