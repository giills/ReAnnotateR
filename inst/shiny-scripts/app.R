# ReAnnotateR Shiny App
# Interactive tool for genomic region annotation and visualization

library(shiny)
library(ReAnnotateR)
library(GenomicRanges)
library(ggplot2)

# UI Definition
ui <- fluidPage(
  titlePanel("ReAnnotateR: Genomic Region Annotation Tool"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Upload"),
      fileInput("bedfile", "Upload BED file",
                accept = c(".bed", ".txt")),
      helpText("Upload a BED file with genomic intervals (tab-separated: chr, start, end)"),
      
      hr(),
      
      h3("Annotation Settings"),
      selectInput("genome", "Reference Genome:",
                  choices = c("hg38" = "hg38", "hg19" = "hg19"),
                  selected = "hg38"),
      
      numericInput("promoter_up", "Promoter upstream (bp):",
                   value = 3000, min = 0, max = 10000, step = 500),
      
      numericInput("promoter_down", "Promoter downstream (bp):",
                   value = 3000, min = 0, max = 10000, step = 500),
      
      hr(),
      
      actionButton("annotate_btn", "Annotate Regions", 
                   class = "btn-primary btn-lg"),
      
      hr(),
      
      h3("Download Results"),
      downloadButton("download_tsv", "Download Annotated Data (TSV)"),
      br(), br(),
      downloadButton("download_bed", "Download Annotated Data (BED)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions",
                 h3("Welcome to ReAnnotateR!"),
                 p("This tool allows you to annotate genomic regions with their feature types 
                   (promoter, exon, intron, intergenic) and identify nearest genes."),
                 
                 h4("Data Format Requirements:"),
                 tags$ul(
                   tags$li("BED format file (tab-separated or space-separated)"),
                   tags$li("Minimum 3 columns: chromosome, start, end"),
                   tags$li("Chromosome names should match reference genome (e.g., chr1, chr2, ...)"),
                   tags$li("0-based or 1-based coordinates accepted")
                 ),
                 
                 h4("Example Data:"),
                 p("Example BED file is available in the package at:"),
                 code("system.file('extdata', 'example.bed', package = 'ReAnnotateR')"),
                 
                 h4("How to Use:"),
                 tags$ol(
                   tags$li("Upload your BED file using the file input"),
                   tags$li("Select reference genome (hg38 or hg19)"),
                   tags$li("Adjust promoter region settings if needed"),
                   tags$li("Click 'Annotate Regions'"),
                   tags$li("Explore results in the tabs"),
                   tags$li("Download annotated data")
                 )
        ),
        
        tabPanel("Data Summary",
                 h3("Uploaded Data"),
                 tableOutput("data_preview"),
                 hr(),
                 verbatimTextOutput("data_summary")
        ),
        
        tabPanel("Feature Composition",
                 h3("Genomic Feature Distribution"),
                 plotOutput("feature_plot", height = "500px"),
                 hr(),
                 tableOutput("feature_table")
        ),
        
        tabPanel("Chromosomal Distribution",
                 h3("Density Across Chromosomes"),
                 plotOutput("density_plot", height = "600px")
        ),
        
        tabPanel("Annotation Summary",
                 h3("Complete Annotation Results"),
                 plotOutput("summary_plot", height = "500px"),
                 hr(),
                 h4("Annotated Regions Table"),
                 DT::dataTableOutput("annotated_table")
        ),
        
        tabPanel("Nearest Genes",
                 h3("Nearest Gene Information"),
                 DT::dataTableOutput("genes_table"),
                 hr(),
                 plotOutput("distance_hist", height = "400px")
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive values to store data
  rv <- reactiveValues(
    raw_gr = NULL,
    annotated_gr = NULL,
    nearest_gr = NULL,
    txdb = NULL
  )
  
  # Load example data on startup
  observe({
    bed_file <- system.file("extdata", "example.bed", package = "ReAnnotateR")
    if (file.exists(bed_file)) {
      rv$raw_gr <- read_bed(bed_file)
    }
  })
  
  # Load BED file
  observeEvent(input$bedfile, {
    req(input$bedfile)
    
    tryCatch({
      rv$raw_gr <- read_bed(input$bedfile$datapath)
      showNotification("BED file loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
    })
  })
  
  # Perform annotation
  observeEvent(input$annotate_btn, {
    req(rv$raw_gr)
    
    withProgress(message = "Annotating regions...", value = 0, {
      tryCatch({
        # Load appropriate TxDb
        incProgress(0.2, detail = "Loading genome annotations...")
        
        if (input$genome == "hg38") {
          if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
            showNotification("TxDb.Hsapiens.UCSC.hg38.knownGene package required!", 
                             type = "error")
            return(NULL)
          }
          rv$txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        } else {
          if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
            showNotification("TxDb.Hsapiens.UCSC.hg19.knownGene package required!", 
                             type = "error")
            return(NULL)
          }
          rv$txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        }
        
        # Annotate regions
        incProgress(0.4, detail = "Annotating features...")
        rv$annotated_gr <- suppressWarnings(
          annotate_regions(rv$raw_gr, rv$txdb,
                           promoter_up = input$promoter_up,
                           promoter_down = input$promoter_down)
        )
        
        # Find nearest genes
        incProgress(0.7, detail = "Finding nearest genes...")
        rv$nearest_gr <- nearest_gene(rv$annotated_gr, rv$txdb)
        
        incProgress(1, detail = "Complete!")
        showNotification("Annotation completed successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Annotation error:", e$message), type = "error")
      })
    })
  })
  
  # Data preview
  output$data_preview <- renderTable({
    req(rv$raw_gr)
    df <- as.data.frame(rv$raw_gr)
    head(df, 10)
  })
  
  # Data summary
  output$data_summary <- renderPrint({
    req(rv$raw_gr)
    cat("Number of regions:", length(rv$raw_gr), "\n")
    cat("Chromosomes:", paste(unique(as.character(seqnames(rv$raw_gr))), collapse = ", "), "\n")
    cat("Total span:", sum(width(rv$raw_gr)), "bp\n")
  })
  
  # Feature composition plot
  output$feature_plot <- renderPlot({
    req(rv$annotated_gr)
    plot_feature_composition(rv$annotated_gr)
  })
  
  # Feature table
  output$feature_table <- renderTable({
    req(rv$annotated_gr)
    as.data.frame(table(mcols(rv$annotated_gr)$feature_type))
  }, colnames = TRUE)
  
  # Chromosomal density plot
  output$density_plot <- renderPlot({
    req(rv$annotated_gr)
    plot_chromosomal_density(rv$annotated_gr)
  })
  
  # Annotation summary plot
  output$summary_plot <- renderPlot({
    req(rv$annotated_gr)
    plot_annotation_summary(rv$annotated_gr)
  })
  
  # Annotated table
  output$annotated_table <- DT::renderDataTable({
    req(rv$annotated_gr)
    df <- as.data.frame(rv$annotated_gr)
    DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE))
  })
  
  # Genes table
  output$genes_table <- DT::renderDataTable({
    req(rv$nearest_gr)
    df <- data.frame(
      chromosome = as.character(seqnames(rv$nearest_gr)),
      start = start(rv$nearest_gr),
      end = end(rv$nearest_gr),
      feature_type = mcols(rv$nearest_gr)$feature_type,
      nearest_gene = mcols(rv$nearest_gr)$nearest_gene,
      distance = mcols(rv$nearest_gr)$distance_to_gene
    )
    DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE))
  })
  
  # Distance histogram
  output$distance_hist <- renderPlot({
    req(rv$nearest_gr)
    distances <- mcols(rv$nearest_gr)$distance_to_gene
    df <- data.frame(distance = distances)
    
    ggplot(df, aes(x = distance)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black") +
      scale_x_log10() +
      labs(title = "Distribution of Distances to Nearest Genes",
           x = "Distance (bp, log scale)",
           y = "Count") +
      theme_minimal()
  })
  
  # Download handlers
  output$download_tsv <- downloadHandler(
    filename = function() {
      paste0("annotated_regions_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(rv$nearest_gr)
      export_results(rv$nearest_gr, file, format = "tsv")
    }
  )
  
  output$download_bed <- downloadHandler(
    filename = function() {
      paste0("annotated_regions_", Sys.Date(), ".bed")
    },
    content = function(file) {
      req(rv$nearest_gr)
      export_results(rv$nearest_gr, file, format = "bed")
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)