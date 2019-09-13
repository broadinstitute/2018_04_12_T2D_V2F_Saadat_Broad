suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))

results_file <- file.path("results", "fat_globule_analysis_results.tsv")
all_results_df <- readr::read_tsv(results_file, col_types=readr::cols())
all_results_df <- all_results_df %>%
    dplyr::mutate(bodipy_color = ifelse(
        grepl("BODIPY", all_results_df$cp_feature),
        'bodipy feature', "not bodipy")
                 )

head(all_results_df, 3)

length(unique(all_results_df$cp_feature))

label_logic <- abs(all_results_df$neg_log_10_p) > 3 | all_results_df$t_stat < -5

append_day <- function(string) paste("Day:", string)
append_ffa <- function(string) paste0("FFA: ", string)

ggplot(all_results_df,
       aes(x=t_stat, y=neg_log_10_p)) +
    geom_point(aes(color = bodipy_color,
                   size = bodipy_color)) +
    xlab("T Stat") +
    ylab("-log10 p value") +
    facet_grid(diff_day~FFA,
               labeller = labeller(diff_day = as_labeller(append_day),
                                   FFA = as_labeller(append_ffa))) +
    scale_size_manual(name = "", 
                      values = c("bodipy feature" = 1, "not bodipy" = 0.02)) +
    theme_bw() +
    geom_text_repel(data = subset(all_results_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.6,
                    point.padding = 0.3,
                    segment.size = 0.2,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    aes(label = cp_feature,
                        x = t_stat,
                        y = neg_log_10_p)) +
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

out_file <- file.path("figures", "volcano_fat_globule_viz.png")
ggsave(out_file, height = 6, width = 7, dpi = 300)

label_logic <- (all_results_df$bodipy_color == "bodipy feature" &
                all_results_df$FFA == 0 &
                all_results_df$t_stat > 5 &
                all_results_df$diff_day == 14) | 
    (all_results_df$bodipy_color == "bodipy feature" &
     all_results_df$FFA == 0 &
     all_results_df$t_stat < -6 &
     all_results_df$diff_day == 14) |
(all_results_df$bodipy_color == "bodipy feature" &
                all_results_df$FFA == 1 &
                all_results_df$t_stat > 5 &
                all_results_df$diff_day == 14) | 
    (all_results_df$bodipy_color == "bodipy feature" &
     all_results_df$FFA == 1 &
     all_results_df$t_stat < -6 &
                all_results_df$diff_day == 14)

ggplot(all_results_df,
       aes(x = diff_day, y = t_stat, group = cp_feature)) +
    geom_point(aes(color = bodipy_color,
                   alpha = bodipy_color,
                   size = bodipy_color)) +
    geom_line(aes(alpha = bodipy_color)) +
    ylab("T Stat") +
    xlab("Differentiation Day") +
    facet_wrap(~FFA,
               labeller = labeller(FFA = as_labeller(append_ffa))) +
    scale_color_discrete(name = "") +
    scale_alpha_manual(name = "",
                       values = c("bodipy feature" = 1,
                                  "not bodipy" = 0.1)) +
    scale_size_manual(name = "", 
                      values = c("bodipy feature" = 1,
                                 "not bodipy" = 0.1)) +
    theme_bw() +
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    geom_text_repel(data = subset(all_results_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.6,
                    point.padding = 0.3,
                    segment.size = 0.2,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    aes(label = cp_feature,
                        x = diff_day,
                        y = t_stat))

out_file <- file.path("figures", "trajectory_fat_globule_viz.png")
ggsave(out_file, height = 5, width = 7, dpi = 300)
