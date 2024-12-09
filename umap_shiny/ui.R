library(shiny)

metadata_choices <- c("Plate", "Well", "Cell_Line", "Patient", "Day", "Batch",
                      "T2D Bin", "T2D Quantile", "HOMA-IR Percentile",
                      "HOMA-IR Rank", "IID", "Category", "None")

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    # Application title
    titlePanel("Type 2 Diabetes UMAP Explorer"),
    
    # Sidebar with interactive layout
    sidebarLayout(

      sidebarPanel(
        helpText("Explore UMAP dimensionality reduction results from Cell Painting data"),
        selectInput("Metadata_Color",
                    label = "Select Metadata to Color",
                    choices = metadata_choices,
                    selected = "Cell_Line"),
        selectInput("Metadata_Shape",
                    label = "Select Metadata to Shape",
                    choices = metadata_choices,
                    selected = "None"),
        selectInput("Differentation_Day",
                    label = "Select Differentiation Day",
                    choices = c("All", "0", "1", "2", "3", "7", "8", "10", "14", "15"),
                    selected = "All"),
        radioButtons(inputId = "filetype",
                     label = "Select file type to download",
                     choices = list("png", "pdf")),
        downloadButton(outputId = "download",
                       label = "Download Plot")
        ),
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("t2d_umap",
                   height = 600,
                   width = 800,
                   click = clickOpts(id = "plot_click"),
                   brush = brushOpts(id = "plot_brush"))
        )
    ),
    fluidRow(
      column(width = 6,
             h4("Points Near Click"),
             verbatimTextOutput("click_info")),
      column(width = 6,
             h4("Click and Drag Points"),
             verbatimTextOutput("brush_info"))
      )
    )
  )
