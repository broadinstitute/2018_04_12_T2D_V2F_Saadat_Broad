save_figure <- function(main_figure,
                        file_base,
                        extensions = c(".png", ".pdf", ".svg"),
                        dpi = 300,
                        height = 6,
                        width = 8) {
  # Save figure given extensions
  #
  # Arguments:
  # main_figure - the cowplot or ggplot object
  # file_base - the name of the file without extensions
  # extensions - character vector of file extensions to save
  # height - height of plot
  # width - width of plot
  #
  # Output:
  # Will save plots to file

  for (extension in extensions) {
    ggsave(main_figure,
           filename = paste0(file_base, extension),
           dpi = dpi,
           height = height,
           width = width)
  }
}
