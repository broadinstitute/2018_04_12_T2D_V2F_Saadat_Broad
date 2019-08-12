
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))

set.seed(123)

batch_id <- "2019_06_11_Batch2"
backend_dir <- file.path("..", "..", "backend", batch_id)

profile_files <- list.files(backend_dir,
                            full.names = TRUE,
                            recursive = TRUE,
                            pattern = "_variable_selected.csv")

profile_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_Plate = readr::col_character(),
  Metadata_Well = readr::col_character(),
  Metadata_Assay_Plate_Barcode = readr::col_character(),
  Metadata_Plate_Map_Name = readr::col_character(),
  Metadata_well_position = readr::col_character(),
  Metadata_cell_line = readr::col_character(),
  Metadata_condition_O2 = readr::col_character()
)

profile_df <- purrr::map_df(
    profile_files,
    readr::read_csv,
    col_types = readr::cols()
)

dim(profile_df)
head(profile_df, 2)

# Output Combined Profiles
file <- file.path("data", paste0("merged_profiles_", batch_id, ".tsv.gz"))
readr::write_tsv(profile_df, file)

# Separate different cell profiler data
cp_features <- colnames(profile_df) %>%
    stringr::str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

length(cp_features)

cp_metadata <- colnames(profile_df) %>%
    stringr::str_subset("^Metadata_")

length(cp_metadata)

# Create a metadata dictionary and dummy variable "group_id"
# "group_id" distinguishes each separate condition including cell line
# "condition_group_id" distinguishes separate conditions ignoring cell line
metadata_df <- profile_df %>%
    dplyr::select(cp_metadata) %>%
    dplyr::mutate(dictionary_id = paste0("id_", dplyr::row_number()),
                  group_id = group_indices(.,
                                           Metadata_cell_line,
                                           Metadata_condition_O2,
                                           Metadata_treatment))

tail(metadata_df)

table(
    metadata_df$Metadata_Plate,
    metadata_df$Metadata_treatment,
    metadata_df$Metadata_condition_O2,
    metadata_df$Metadata_cell_line
)

# Create a dataframe of variables for each group
group_id_df <- metadata_df %>%
    dplyr::select(
        group_id,
        Metadata_cell_line,
        Metadata_condition_O2,
        Metadata_treatment
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(group_id)

dim(group_id_df)
head(group_id_df)

cor_df <- profile_df %>%
    dplyr::select(cp_features) %>%
    t() %>%
    cor() %>%
    dplyr::as_tibble() %>%
    magrittr::set_colnames(metadata_df$dictionary_id)

cor_melt_df <- metadata_df %>%
    dplyr::select(-group_id) %>%
    dplyr::bind_cols(
        replace(cor_df,
                lower.tri(cor_df, TRUE), NA)
    ) %>%
    dplyr::select(-cp_metadata) %>%
    reshape2::melt(id.vars = 'dictionary_id',
                   variable.name = 'correlation_id', 
                   value.name = "pearson_cor",
                   na.rm = TRUE) %>%
    tibble::remove_rownames()

dim(cor_melt_df)
head(cor_melt_df)

# Map group IDs and condition IDs onto the correlation dataframe
# We are interested in correlations between specific groups
cor_group_df <- cor_melt_df %>%
    dplyr::inner_join(
        metadata_df %>%
        select(dictionary_id,
               group_id),
        by = 'dictionary_id'
    ) %>%
    dplyr::rename(pair_a = group_id) %>%
    dplyr::inner_join(
        metadata_df %>%
        select(dictionary_id,
               group_id),
        by = c('correlation_id' = 'dictionary_id')
    ) %>%
    dplyr::rename(pair_b = group_id,
                  pair_a_id = dictionary_id,
                  pair_b_id = correlation_id)

dim(cor_group_df)
head(cor_group_df)

# Remove self correlations and determine median correlation between all groups
# Also create a variable that represents correlations across cell lines within
# the same condition. This variable will be used as the null distribution.
cor_group_df <- cor_group_df %>%
    dplyr::mutate(
        within_group_cor =
            as.numeric(cor_group_df$pair_a == cor_group_df$pair_b),
    ) %>%
    dplyr::filter(cor_group_df$pair_a_id != cor_group_df$pair_b_id) %>%
    dplyr::group_by(
        pair_a,
        pair_b
    ) %>%
    dplyr::mutate(median_cor = median(pearson_cor)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(median_cor))

dim(cor_group_df)
head(cor_group_df)

# Join Replicate Correlations and Null Distribution Correlations
within_group_cor_df <- cor_group_df %>%
    dplyr::filter(within_group_cor == 1) %>%
    dplyr::group_by(pair_b) %>%
    dplyr::mutate(pair_b_median_cor = median(pearson_cor),
                  null_data = "Replicate Correlation") %>%
    dplyr::arrange(desc(pair_b_median_cor)) %>%
    dplyr::ungroup()

null_group_cor_df <- cor_group_df %>%
    dplyr::filter(within_group_cor == 0) %>%
    dplyr::group_by(pair_b) %>%
    dplyr::mutate(pair_b_median_cor = median(pearson_cor),
                  null_data = "Non Replicate Correlation") %>%
    dplyr::arrange(desc(pair_b_median_cor)) %>%
    dplyr::ungroup()

full_plot_ready <- within_group_cor_df %>%
    dplyr::bind_rows(
        null_group_cor_df
    )

dim(full_plot_ready)
head(full_plot_ready)

# Merge plot ready data with info on group ID
full_plot_ready <- full_plot_ready %>%
    dplyr::left_join(group_id_df, by = c("pair_b" = "group_id"))

head(full_plot_ready, 2)

# Update Factors for plotting
full_plot_ready$Metadata_condition_O2 <-
    factor(full_plot_ready$Metadata_condition_O2,
           levels = c("4", "21"))

append_02 <- function(string) paste("O2:", string)

cor_gg <- ggplot(full_plot_ready,
                 aes(x = pearson_cor,
                     fill = null_data)) +
    geom_density(alpha = 0.8) +
    ylab("Patient") +
    xlab("Profile Correlation (Pearson)") +
    theme_bw() +
    facet_grid(Metadata_treatment~Metadata_condition_O2,
               scales="free_y",
               labeller = labeller(Metadata_condition_O2 = as_labeller(append_02))) +
    scale_fill_manual(name = "",
                      values = c("#FFC107", "#004D40")) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 7),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

print(cor_gg)

file_base <- file.path("figures", paste0("replicate_correlation_batch_", batch_id))
for (extension in c('.png', '.pdf')) {
    ggsave(cor_gg,
           filename = paste0(file_base, extension),
           height = 10,
           width = 8)
}

# Perform KS tests between real and null distributions
all_results <- list()
for (group_id in unique(full_plot_ready$pair_b)) {
    full_plot_group_df = full_plot_ready %>%
        dplyr::filter(pair_b == group_id)

    replicate_corr_df <- full_plot_group_df %>%
        dplyr::filter(within_group_cor == 1)
    null_corr_df <- full_plot_group_df %>%
        dplyr::filter(within_group_cor == 0)

    ks_result = ks.test(x = replicate_corr_df$pearson_cor,
                        y = null_corr_df$pearson_cor,
                        alternative = "less")

    k_stat = as.numeric(ks_result$statistic)
    k_p = as.numeric(ks_result$p.value)
    all_results[[group_id]] <- c(group_id, k_stat, k_p, -log10(k_p), nrow(replicate_corr_df))
}

ks_result_df <- dplyr::as_tibble(do.call(rbind, all_results))
colnames(ks_result_df) <- c("group_id", "ks_stat", "ks_p_value", "ks_log_10_p", "num_replicates")

ks_result_df <- ks_result_df %>% dplyr::arrange(desc(as.numeric(paste(ks_log_10_p))))

dim(ks_result_df)
head(ks_result_df)

ks_result_df$group_id <- as.numeric(ks_result_df$group_id)

final_results_df <- full_plot_ready %>%
    dplyr::left_join(ks_result_df,
                     by = c("pair_b" = "group_id")) %>%
    dplyr::group_by(
        Metadata_cell_line,
        Metadata_condition_O2,
        Metadata_treatment
    ) %>%
    dplyr::mutate(
        median_cell_ks = median(ks_stat)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
        Metadata_cell_line,
        Metadata_condition_O2,
        Metadata_treatment,
        ks_stat,
        ks_p_value
    ) %>%
    dplyr::distinct()

head(final_results_df, 2)

ks_test_gg <- ggplot(final_results_df) +
    geom_bar(aes(y = ks_stat,
                 x = as.factor(Metadata_treatment),
                 fill = Metadata_condition_O2),
                stat = "identity",
            position = "dodge") +
    xlab("Treatment") +
    ylab("KS Statistic") +
    scale_fill_manual(name = "O2 Condition",
                       values = c("4" = "#a6cee3", "21" = "#d95f02"),
                       labels = c("4" = "4", "21" = "21")) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 8,
                                     angle = 45,
                                     vjust = 0.5),
          axis.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 9),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

ks_test_gg

file_base <- file.path("figures", paste0("replicate_correlation_kstest", "_", batch_id))
for (extension in c('.png', '.pdf')) {
    ggsave(ks_test_gg,
           filename = paste0(file_base, extension),
           height = 4,
           width = 6)
}
