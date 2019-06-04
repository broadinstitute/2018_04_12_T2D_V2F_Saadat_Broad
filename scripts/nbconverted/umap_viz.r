
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(umap))

set.seed(12345)

file <- file.path("data", "combined_normalized_variable_selected.tsv")

cp_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_Plate = readr::col_character(),
  Metadata_Well = readr::col_character(),
  Metadata_Assay_Plate_Barcode = readr::col_character(),
  Metadata_Plate_Map_Name = readr::col_character(),
  Metadata_well_position = readr::col_character(),
  Metadata_cell_line = readr::col_character(),
  Metadata_patient = readr::col_character(),
  Metadata_FFA = readr::col_character(),
  Metadata_diff_day = readr::col_character()
)

df <- readr::read_tsv(file, col_types = cp_cols)

dim(df)
head(df, 2)

metadata_df <- df %>%
    dplyr::select(dplyr::starts_with("Metadata_"))

cp_df <- df %>%
    dplyr::select(-dplyr::starts_with("Metadata_"))

# Apply UMAP
cp_umap <- umap(as.matrix(cp_df))
cp_umap_df <- cp_umap$layout %>%
    dplyr::as_tibble()

colnames(cp_umap_df) <- c("umap_x", "umap_y")

# Merge with metadata
cp_umap_df <- cp_umap_df %>%
    dplyr::bind_cols(metadata_df)

head(cp_umap_df, 2)

# Write umap output
file <- file.path("results", "umap_with_metadata.tsv")
readr::write_tsv(cp_umap_df, file)

# Reorder some factors
cp_umap_df$Metadata_diff_day <-
    factor(cp_umap_df$Metadata_diff_day, levels = c("0", "1", "2", "3", "7", "10", "14", "15", "15+iso"))
cp_umap_df$Metadata_patient <-
    factor(cp_umap_df$Metadata_patient,
           levels = c("PAC_164", "PAC_246", "PAC_261", "PAC_266", "hBAT", "hWAT", "SGBS"))

# Set Constants
legend_text_size = 7
legend_title_size = 8
legend_key_height = 0.18
axis_title_size = 8
axis_text_size = 7
geom_point_size = 0.8

plate_gg <-
    ggplot(cp_umap_df, aes(umap_x, umap_y)) +
        geom_point(aes(fill = Metadata_Plate),
                   size = geom_point_size,
                   alpha = 0.5,
                   color = "black",
                   pch = 21) +
        theme_bw() +
        xlab("UMAP (x)") +
        ylab("UMAP (y)") +
        scale_fill_viridis_d(name = "Plate", option = "plasma") +
        guides(fill = guide_legend(keyheight = legend_key_height,
                                   default.unit = "inch",
                                   override.aes = list(alpha = 1,
                                                       size = 1))) +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.title = element_text(size = axis_title_size),
              axis.text = element_text(size = axis_text_size))

plate_gg

patient_gg <-
    ggplot(cp_umap_df, aes(umap_x, umap_y)) +
        geom_point(aes(fill = Metadata_patient),
                   size = geom_point_size,
                   alpha = 0.5,
                   pch = 21,
                   color = "black") +
        theme_bw() +
        xlab("UMAP (x)") +
        ylab("UMAP (y)") +
        scale_fill_manual(name = "Patient",
                           values = c("#66c2a5",
                                      "#fc8d62",
                                      "#8da0cb",
                                      "#e78ac3",
                                      "#a6d854",
                                      "#ffd92f",
                                      "#e5c494")) +
        guides(fill = guide_legend(keyheight = legend_key_height,
                                   default.unit = "inch",
                                   override.aes = list(alpha = 1,
                                                       size = 1))) +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.title = element_text(size = axis_title_size),
              axis.text = element_text(size = axis_text_size))

patient_gg

day_gg <- 
    ggplot(cp_umap_df, aes(umap_x, umap_y)) +
        geom_point(aes(fill = Metadata_diff_day),
                   size = geom_point_size,
                   alpha = 0.5,
                   pch = 21,
                   color = "black") +
        xlab("UMAP (x)") +
        ylab("UMAP (y)") +
        theme_bw() +
        scale_fill_viridis_d(name = "Day") +
        guides(fill = guide_legend(keyheight = legend_key_height,
                                   default.unit = "inch",
                                   override.aes = list(alpha = 1,
                                                       size = 1))) +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.title = element_text(size = axis_title_size),
              axis.text = element_text(size = axis_text_size))

day_gg

ffa_gg <-
    ggplot(cp_umap_df, aes(umap_x, umap_y)) +
        geom_point(aes(fill = Metadata_FFA),
                   size = geom_point_size,
                   alpha = 0.5,
                   pch = 21,
                   color = 'black') +
        xlab("UMAP (x)") +
        ylab("UMAP (y)") +
        theme_bw() +
        scale_fill_manual(name = "Free Fatty Acids",
                          values = c("0" = "#f2f4f7", "1" = "#0a4cb5"),
                          labels = c("0" = "None", "1" = "Added")) +
        guides(fill = guide_legend(keyheight = legend_key_height,
                                   default.unit = "inch",
                                   override.aes = list(alpha = 1,
                                                       size = 1))) +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.title = element_text(size = axis_title_size),
              axis.text = element_text(size = axis_text_size))

ffa_gg

main_plot <- (
    cowplot::plot_grid(
        plate_gg,
        patient_gg,
        day_gg,
        ffa_gg,
        labels = c("a", "b", "c", "d"),
        ncol = 2,
        nrow = 2,
        align = "v"
    )
)

main_plot

for(extension in c('.png', '.pdf')) {
    sup_file <- paste0("umap_metadata", extension)
    sup_file <- file.path("figures", sup_file)
    cowplot::save_plot(filename = sup_file,
                       plot = main_plot,
                       base_height = 130,
                       base_width = 200,
                       unit = "mm")
}
