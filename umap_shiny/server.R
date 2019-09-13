suppressMessages(library(shiny))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# Load UMAP results
umap_file <- file.path("data", "combined_batch1_batch3_umap_with_metadata.tsv")

umap_cols <- readr::cols(
  .default = readr::col_double(),
  Batch = readr::col_character(),
  Plate = readr::col_character(),
  Well = readr::col_character(),
  Cell_Line = readr::col_character(),
  Patient = readr::col_character(),
  IID = readr::col_character(),
  `T2D Bin`	= readr::col_character(),
  Category = readr::col_character(),
  Day = readr::col_character()
)

umap_df <- readr::read_tsv(umap_file, col_types = umap_cols)

shinyServer(function(input, output) {
  
  output$plot_click <- renderPrint({
    cat("input$plot_click:\n")
    str(input$Metadata_Color)
  })
  
  metadata_color_ <- reactive({
    paste(input$Metadata_Color)
  })
  
  metadata_shape_ <- reactive({
    paste(input$Metadata_Shape)
  })
  
  diff_day_select_ <- reactive({
    paste(input$Differentation_Day)
  })

  output$t2d_umap <- renderPlot(
    {
      metadata_color <- metadata_color_()
      metadata_shape <- metadata_shape_()
      diff_day_subset <- diff_day_select_()
      
      # Determine which data to plot
      if (diff_day_subset == "All") {
        plot_df <- umap_df
      } else {
        plot_df <- umap_df %>% dplyr::filter(Day == diff_day_subset)
      }
      
      p <- ggplot(plot_df, aes(x, y)) + 
        xlab("UMAP (x)") +
        ylab("UMAP (y)") +
        theme_bw()
      
      if ((metadata_color != "None") & (metadata_shape == "None")) {
        p <- p + geom_point(aes_string(color = metadata_color),
                            size = 2, alpha = 0.6)
      } else if ((metadata_color != "None") & (metadata_shape != "None")) {
        p <- p + geom_point(aes_string(color = metadata_color,
                                       shape = metadata_shape),
                            size = 2, alpha = 0.6)
      } else if ((metadata_color == "None") & (metadata_shape != "None")) {
        p <- p + geom_point(aes_string(shape = metadata_shape),
                            size = 2, alpha = 0.6)
      } else {
        p <- p + geom_point(size = 2, alpha = 0.6)
      }

      p
    })
  
  output$click_info <- renderPrint({
    nearPoints(umap_df, input$plot_click, addDist = FALSE)
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(umap_df, input$plot_brush)
  })
  
  output$download <- downloadHandler(
    filename = function() {
      base <- paste0("umap_color_", metadata_color_(), "_shape_", metadata_shape_())
      paste0(base, ".", input$filetype)
    },
    
    content = function(file) {
      # same as renderplot
      metadata_color <- metadata_color_()
      metadata_shape <- metadata_shape_()
      diff_day_subset <- diff_day_select_()
      
      # Determine which data to plot
      if (diff_day_subset == "All") {
        plot_df <- umap_df
      } else {
        plot_df <- umap_df %>% dplyr::filter(Day == diff_day_subset)
      }
      
      p <- ggplot(plot_df, aes(x, y)) + 
        xlab("UMAP (x)") +
        ylab("UMAP (y)") +
        theme_bw()
      
      if ((metadata_color != "None") & (metadata_shape == "None")) {
        p <- p + geom_point(aes_string(color = metadata_color),
                            size = 1.2, alpha = 0.6)
      } else if ((metadata_color != "None") & (metadata_shape != "None")) {
        p <- p + geom_point(aes_string(color = metadata_color,
                                       shape = metadata_shape),
                            size = 1.2, alpha = 0.6)
      } else if ((metadata_color == "None") & (metadata_shape != "None")) {
        p <- p + geom_point(aes_string(shape = metadata_shape),
                            size = 1.2, alpha = 0.6)
      } else {
        p <- p + geom_point(size = 1.2, alpha = 0.6)
      }
      p + ggsave(file, width = 6, height = 5)
    }
  )
})