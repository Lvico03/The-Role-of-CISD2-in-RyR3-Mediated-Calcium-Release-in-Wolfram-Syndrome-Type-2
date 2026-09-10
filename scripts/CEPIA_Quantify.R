
library(tidyr)
library(purrr)
library(data.table)
library(plotly)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(plyr)
library(dplyr)
library(car)
library(readxl)
library(janitor)
library(data.table)
library(ggsci)
library(doBy)
library(Rmisc)

files <- list.files(pattern = "\\.csv$", full.names=TRUE)
transposed_list <- lapply(files, function(f) {
  df <- read.csv(f, row.names = 1)
  t_df <- t(df)
  as.data.frame(t_df)
})
common_cols <- Reduce(intersect, lapply(transposed_list, colnames))
transposed_list <- lapply(transposed_list,function(df) df[, common_cols, drop = FALSE])
merged_df <- bind_rows(transposed_list)
bind_rows(transposed_list)
write.csv(merged_df, "20250724_Hek_RyR3_CISD2_KO_QuantifyResponse.csv", row.names=FALSE)

#Quantify Response
library(tidyverse)

files <- list.files(pattern = "\\.csv$")
df_all <- files %>%
  map_dfr(read_csv) %>%
  mutate(
    Day = str_extract(ID, "^\\d{8}"),                  
    ID_short = str_remove(ID, "^\\d{8}_"),             
    Day = factor(Day, sort(unique(Day)),
                 labels = paste0("Day", seq_along(sort(unique(Day)))))     
  )

df_summary <- df_all %>%
  group_by(ID_short, Genotype, Day) %>%
  summarise(
    Max_amplitude = mean(Max_amplitude, na.rm = TRUE),
    AUC = mean(AUC, na.rm = TRUE),
    t_MAX = mean(t_MAX, na.rm = TRUE),
    Basal = mean(Basal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Genotype = factor (Genotype, levels = c("RyR3_WT", "RyR3_CISD2_KO")))


#Anova for AUC after summarize
fit_summary <- aov(AUC ~ Genotype, data = df_summary)
cooksd_summary <- cooks.distance(fit_summary)
sample_size_summary <- nrow(df_summary)

#Cook's distance
plot(cooksd_summary, pch = "*", cex = 2,
     main = "Influential Obs by Cook's distance AUC summarized")
abline( h = 4/ sample_size_summary, col = "red")
text(x = 1:length(cooksd_summary), 
     y = cooksd_summary,
     labels = ifelse(cooksd_summary > 4 / sample_size_summary, 
                     df_summary$ID_short, ""), 
     col = "red", pos = 4, cex = 0.8)

#Outlier removal

influential_summary <- which(cooksd_summary > 4 / sample_size_summary)
df_summary_clean <- df_summary[-influential_summary, ]


# Norm
wt_means <- df_summary %>%
  filter(Genotype == "RyR3_WT") %>%
  select(Day, Max_amplitude, AUC, Basal, t_MAX) %>%
  group_by(Day) %>%
  summarise(
    Max_amplitude_WT = mean(Max_amplitude, na.rm = TRUE),
    AUC_WT = mean(AUC, na.rm = TRUE),
    Basal_WT = mean(Basal, na.rm = TRUE),
    t_MAX_WT = mean(t_MAX, na.rm = TRUE),
    .groups = "drop"
  )

df_norm <- df_summary %>%
  left_join(wt_means, by = "Day") %>%
  mutate(
    Max_amplitude_norm = Max_amplitude / Max_amplitude_WT,
    AUC_norm = AUC / AUC_WT,
    Basal_norm = Basal / Basal_WT,
    t_MAX_norm = t_MAX / t_MAX_WT
  )


#BAR PLOT
make_plot <- function(data, variable, y_label) {
  ggplot(data, aes(x = Genotype, y = .data[[variable]], fill = Genotype)) +
    stat_summary(fun = mean, geom = "bar", color = "black", width = 0.3,
                 position = position_dodge(width = 0.8)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.05,
                 position = position_dodge(width = 0.8)) +
    geom_point(aes(shape = Day), size = 3,
               position = position_jitterdodge(jitter.width = 0, dodge.width = 0.4)) +
    labs(title = paste(y_label, "by Genotype"),
         y = y_label, x = "") +
    theme_minimal() +
    scale_shape_manual(values = setNames(16:(16 + nlevels(data$Day) - 1), levels(data$Day))) +
    scale_fill_manual(values = c("RyR3_WT" = "grey", "RyR3_CISD2_KO" = "green4")) +
    theme(legend.position = "right")
}

# Not normalized Bar Plots
p1 <- make_plot(df_summary, "Max_amplitude", "Max Amplitude")
print(p1)

p2 <- make_plot(df_summary, "AUC", "AUC")
print(p2)

p3 <- make_plot(df_summary, "Basal", "Basal")
print(p3)

p4 <- make_plot(df_summary, "t_MAX", "t_MAX")
print(p4)

# Normalized Bar Plots
p1_norm <- make_plot(df_norm, "Max_amplitude_norm", "Normalized Max Amplitude")
print(p1_norm)

p2_norm <- make_plot(df_norm, "AUC_norm", "Normalized AUC")
print(p2_norm)

p3_norm <- make_plot(df_norm, "Basal_norm", "Normalized Basal")
print(p3_norm)

p4_norm <- make_plot(df_norm, "t_MAX_norm", "Normalized t_MAX")
print(p4_norm)

#VIOLIN PLOT

plot_violin = ggplot(df_norm, aes(x=Genotype, y=AUC_norm, group = Genotype)) +
  geom_violin(data = df_norm, aes(fill=Genotype, x=Genotype, y=AUC_norm), 
              scale="area", alpha=0.2) +
  geom_sina(data = df_filter, aes(colour = Genotype, x=Genotype, y=AUC), 
            alpha=0.3, size=1.5) +  
  geom_point(aes(shape = Day), 
             color = "black", alpha = 0.9, size = 2) + 
  stat_summary(fun.data = median_hilow, geom = "errorbar", 
               width = 0.2, color = "black", size = 0.6) +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.1, fatten = 0, size = 0.6, color = "black") +
  scale_colour_manual(values = c("green4", "grey")) +
  scale_fill_manual(values = c("green4", "grey")) +
  xlab("") +
  ylab("Basal calcium") +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 13), 
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 14, face = "bold"), 
    strip.text.x = element_text(size = 12, color = "black", face = "bold")
  )


plot_violin