suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(umap))

set.seed(12345)

file <- file.path("data", "batch1_batch3_combined_normalized_variable_selected.tsv")

cp_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_Batch = readr::col_character(),
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

cp_umap_df <- cp_umap_df %>%
    dplyr::select(umap_x, umap_y, Metadata_Plate, Metadata_Well, Metadata_cell_line,
                  Metadata_patient, Metadata_FFA, Metadata_diff_day, Metadata_Batch) %>%
    dplyr::rename(x = umap_x, y = umap_y, Plate = Metadata_Plate,
                  Well = Metadata_Well, Cell_Line = Metadata_cell_line,
                  Patient = Metadata_patient, FFA = Metadata_FFA,
                  Day = Metadata_diff_day, Batch = Metadata_Batch)

head(cp_umap_df, 2)

# Write umap output
file <- file.path("umap_shiny", "data", "combined_batch1_batch3_umap_with_metadata.tsv")
readr::write_tsv(cp_umap_df, file)

patient_gg <- ggplot(cp_umap_df, aes(x, y)) +
    geom_point(aes(fill = Patient),
               size = 1.2,
               alpha = 0.5,
               color = "black",
               pch = 21) +
    theme_bw() +
    scale_fill_discrete(name = "Patient") +
    xlab("UMAP (x)") +
    ylab("UMAP (y)")

patient_gg

batch_gg <- ggplot(cp_umap_df, aes(x, y)) +
    geom_point(aes(fill = Batch),
               size = 1.2,
               alpha = 0.5,
               color = "black",
               pch = 21) +
    theme_bw() +
    scale_fill_discrete(name = "Batch",
                        labels = c("batch_one" = "1",
                                   "batch_three" = "3")) +
    xlab("UMAP (x)") +
    ylab("UMAP (y)")
 
batch_gg

day_gg <- ggplot(cp_umap_df, aes(x, y)) +
    geom_point(aes(fill = Day),
               size = 1.2,
               alpha = 0.5,
               color = "black",
               pch = 21) +
    theme_bw() +
    scale_fill_discrete(name = "Cell Line") +
    xlab("UMAP (x)") +
    ylab("UMAP (y)")
 
day_gg

ffa_gg <- ggplot(cp_umap_df, aes(x, y)) +
    geom_point(aes(fill = FFA),
               size = 1.2,
               alpha = 0.5,
               color = "black",
               pch = 21) +
    theme_bw() +
    scale_fill_discrete(name = "FFA") +
    xlab("UMAP (x)") +
    ylab("UMAP (y)")
 
ffa_gg

main_plot <- (
    cowplot::plot_grid(
        batch_gg,
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
    sup_file <- paste0("umap_metadata_batch1_batch3_combined", extension)
    sup_file <- file.path("figures", sup_file)
    cowplot::save_plot(filename = sup_file,
                       plot = main_plot,
                       base_height = 130,
                       base_width = 200,
                       unit = "mm")
}

cp_umap_df$Day <- dplyr::recode(cp_umap_df$Day, "15+iso" = "15")
cp_umap_df$Day <- factor(cp_umap_df$Day, levels = sort(as.numeric(paste(unique(cp_umap_df$Day)))))

ggplot(cp_umap_df, aes(x, y)) +
    geom_point(aes(color = Cell_Line,
                   size = as.numeric(paste(Day)),
                   shape = Batch),
               alpha = 0.3) +
    theme_bw() +
    scale_size_continuous(name = "Day", range = c(0.5, 2.5)) +
    scale_color_discrete(name = "Cell Line") +
    scale_shape_manual(name = "Batch",
                       values = c(19, 17),
                       labels = c("batch_one" = "1", "batch_three" = "3")) +
    xlab("UMAP (x)") +
    ylab("UMAP (y)") +
    theme(strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 7),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "umap_batch1_batch3_day_line_batch.png")
ggsave(output_file, height = 5, width = 6, dpi = 300)
