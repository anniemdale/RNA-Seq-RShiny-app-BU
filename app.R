library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)

# source the helpers functions to be used here 
source("helpers.R")

#create interface with 4 main tabs with varying sub tabs
ui <- fluidPage(
  mainPanel(
    tabsetPanel(
      tabPanel("Samples", # Main Tab 1 - Sample Information Exploration
               sidebarLayout(
                 sidebarPanel(
                   # Load Sample Data CSV File in sidebar
                   fileInput("datafile", "Choose Sample Info CSV file",
                             accept = c("test/csv", "text/comma-separated-files, .csv")),
                 ),
                 mainPanel(
                   # Main Panel Samples has 3 sub tabs
                   tabsetPanel(
                     # Tab 1 - Summary of Sample Information in a table
                     tabPanel("Summary",tableOutput("sample_summary_table")),
                     # Tab 2 - Sample Information in a sort-able table
                     tabPanel("Table", DTOutput("sample_table")),
                     # Tab 3 - Bar Plots of counts of values in each column
                     tabPanel("Plots",
                              # Side bar to select which column you want to visualize the counts of
                              sidebarLayout(
                                sidebarPanel(
                                  # Drop down menu to select column
                                  selectInput("histogram_column", "Choose column", c("sample","treatment","lifestage","timepoint","period","sex","replicate"), "lifestage")
                                )
                              ,
                              # Main Panel to show the plot
                              mainPanel("Plot", plotOutput("histogram")))
                   )
                 )
               ))),
      tabPanel("Counts", # Main Tab 2 - DESeq Normalized Counts Data Exploration
               sidebarLayout(
                 sidebarPanel(
                   # Load Normalized Counts CSV in the sidebar
                   fileInput("countsfile", "Choose Normalized Counts CSV file",
                             accept = c("test/csv", "text/comma-separated-files, .csv")),
                   # Load Sample Information CSV in the sidebar
                   fileInput("samplefile", "Choose Sample Info CSV file",
                             accept = c("test/csv", "text/comma-separated-files, .csv")),
                   # Sliders to filter genes based on in the sidebar
                   sliderInput("var_slider",
                               "Select Percentile of Variance:",
                               min=0, max=100, value=70),
                   sliderInput("zero_slider",
                               "Select Number of Non-zero Genes:",
                               min=0, max=17696, value=50)
                 ),
                 mainPanel(
                   # Main Panel has 4 sub tabs which update according to the sidebar input
                   tabsetPanel(
                     # Tab 1 - Summary of number and percent of genes passing/failing the filter in a table
                     tabPanel("Summary", tableOutput("counts_summary_table")),
                     # Tab 2 - Scatter plots showing the variance in the genes and number of zero-reads in the genes v. median
                     tabPanel("Scatter", plotOutput("variance_scatter"), plotOutput("zero_scatter")),
                     # Tab 3 - Clustered Heatmap showing the gene counts/reads for each sample
                     tabPanel("Heatmap", plotOutput("heatmap")),
                     # Tab 4 - PCA plot with ability to control PCs plotted against each other
                     tabPanel("PCA",
                              sidebarLayout(
                                sidebarPanel(
                                  # Sliders to select the PCs
                                  sliderInput("first_pc",
                                              "Select PC for x-axis:",
                                              min=1, max=20, value=1),
                                  sliderInput("second_pc",
                                              "Select PC for y-axis:",
                                              min=1, max=20, value=2)
                                ),
                                # Main Panel shows the updated PCA plot based on the selected PCs
                                mainPanel("PCA", plotOutput("pca"))
                              ))
                   )
                 )
               )),
      tabPanel("DE", # Main Tab 3 - DESeq Differential Expression Analysis Results
               sidebarLayout(
                 sidebarPanel(
                   # Upload DESeq Results object as a csv file in the sidebar
                   fileInput("defile", "Choose DESeq2 Results CSV file",
                             accept = c("test/csv", "text/comma-separated-files, .csv")),
                 ),
                 # Main Panel has 2 sub tabs that do not update
                 mainPanel(
                   tabsetPanel(
                     # Tab 1 - Sortable table of the DESeq DE Results
                     tabPanel("DE Results",DTOutput("de_table")),
                     # Tab 2 - Volcano plot showing genes that are DE based on pval < 0.01 & |log2FC| > 0.58
                     tabPanel("Volcano", plotOutput("volcano"),
                     )
                   )
                 ))),
      tabPanel("GSEA", # Main Tab 4 - FGSEA Results and Exploration of Biological Pathways
               sidebarLayout(
                 sidebarPanel(
                   # Upload FGSEA Results File as a csv in the sidebar
                   fileInput("fgseafile", "Choose FGSEA Results CSV file",
                             accept = c("test/csv", "text/comma-separated-files, .csv")),
                   # Slider to adjust pathways being shown based on padj value
                   sliderInput("padj_slider", "Select adjusted p-value:", min=0, max=1, value=0.5)
                 ),
                 # Main Panel has 3 sub tabs that update based on the padj slider
                 mainPanel(
                   tabsetPanel(
                     # Tab 1 - Bar Plot of Top Up/Down Enriched Biological Pathways based on Gene Sets
                     tabPanel("Top Pathways", plotOutput("bar")),
                     # Tab 2 - Downloadable table of either all/positive/negative enriched pathways
                     tabPanel("Table",
                              sidebarLayout(
                                sidebarPanel(
                                  radioButtons("value_button", "Choose NES Pathways:", c("All", "Positive", "Negative"), "All"),
                                  actionButton("loadData", "Load Data"),
                                  downloadButton('download',"Download"),
                                ),
                                mainPanel("Pathways Table", DTOutput("fgsea_table")))),
                     # Tab 3 - Scatter Plot of all pathways based on passing/failing the filter
                     tabPanel("Scatter", plotOutput("fgsea_scatter"))
                     )
                   )
                 )))
    )
  )

# server specifies back end parameters
server <- function(input, output, session) {
  
  # Load sample data
  load_sample_data <- reactive({
    req(input$datafile)              
    data <- read.csv(req(input$datafile$datapath))
    return(data)
  })
  
  # Load normalized counts data
  load_counts_data <- reactive ({
    req(input$countsfile)
    data <- read.csv(req(input$countsfile$datapath))
    return(data)
  })
  
  # Load Sample data again
  load_sample <- reactive ({
    req(input$samplefile)
    data <- read.csv(req(input$samplefile$datapath))
    return(data)
  })
  
  # Load DESeq Reults object data
  load_de_results <- reactive ({
    req(input$defile)
    data <- read.csv(req(input$defile$datapath))
    return(data)
  })
  
  # Load FGSEA Results data
  load_fgsea_results <- reactive ({
    req(input$fgseafile)
    data <- read.csv(req(input$fgseafile$datapath))
    return(data)
  })
  
  # Sample exploration
  # output histogram, calls helper function
  output$histogram <-({
    renderPlot(create_sample_histogram(load_sample_data(), input$histogram_column))
  })
  # output sample summary table
  output$sample_summary_table <- ({
    renderTable(sample_summary(load_sample_data()))
  })
  # output sample table, calls helper function
  output$sample_table <- ({
    renderDT(sample_datatable(load_sample_data()))
  })
  
  # Counts exploration
  # output summary of the normalized counts, calls helper function
  output$counts_summary_table <- ({
    renderTable(count_summary_table(load_counts_data(), input$var_slider, input$zero_slider))
  })
  # output of scatter plot involving variance, calls helper function
  output$variance_scatter <- ({
    renderPlot(count_variance_scatter(load_counts_data(), input$var_slider, input$zero_slider))
  })
  # output of scatter plot involving counting zeros, calls helper function
  output$zero_scatter <- ({
    renderPlot(count_zero_scatter(load_counts_data(), input$var_slider, input$zero_slider))
  })
  # output of clustered heatmap, calls helper function
  output$heatmap <- ({
    renderPlot(count_heatmap(load_counts_data(), input$var_slider, input$zero_slider))
  })
  # output of pca, calls helper function
  output$pca <- ({
    renderPlot(plot_pca(load_counts_data(), load_sample(), input$first_pc, input$second_pc, input$var_slider, input$zero_slider))
  })
  
  # DEG exploration
  # output of DE results table, calls helper function
  output$de_table <- ({
    renderDT(de_datatable(load_de_results()))
  })
  # output of volcano plot to show DE genes, calls helper function
  output$volcano <- ({
    renderPlot(plot_volcano(load_de_results()))
  })
  
  #FGSEA exploration
  # output of top pathways bar plot, calls helper function
  output$bar <- ({
    renderPlot(top_pathways_bar(load_fgsea_results(), input$padj_slider))
  })
  # makes the filtered data based on all/pos/neg selection reactive, calls helper function
  filtered_data <- reactive({
    req(input$loadData)
    isolate({
      fgsea_table(load_fgsea_results(), input$padj_slider, input$value_button)
    })
  })
  # output of fgsea results, using above function
  output$fgsea_table <- renderDT({
    filtered_data()
  })
  # output of downloaded results based on selection, based on above 2 functions
  output$download <- downloadHandler(
    filename = function() {
      paste(input$value_button, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
  # output of fgsea scatter plot, calls helper function
  output$fgsea_scatter <- ({
    renderPlot(fgsea_scatter(load_fgsea_results(), input$padj_slider))
  })
    
}

# allows us to upload larger files
options(shiny.maxRequestSize = 30 * 1024^2)
# call the app
shinyApp(ui = ui, server = server)
