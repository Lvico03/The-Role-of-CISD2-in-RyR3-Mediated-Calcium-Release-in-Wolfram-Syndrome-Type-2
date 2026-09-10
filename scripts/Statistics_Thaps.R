
library(tidyverse)
library(car)

files <- list.files(pattern = "\\.csv$")
df_all <- files %>%
  map_dfr(read_csv) %>%
  mutate(
    Day = str_extract(ID, "^\\d{8}"),                  
    ID_short = str_remove(ID, "^\\d{8}_"),             
    Day = factor(Day, sort(unique(Day)),
                 labels = paste0("Day", seq_along(sort(unique(Day)))))     
  )

#Filtering Day 2, CISD2 KO 1
df_filter <- df_all %>% filter(ID != "20250617_CISD2_KO_1_1")

df_summary <- df_filter %>%
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
wt_means <- df_filter %>%
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

df_norm <- df_filter %>%
  left_join(wt_means, by = "Day") %>%
  mutate(
    Max_amplitude_norm = Max_amplitude / Max_amplitude_WT,
    AUC_norm = AUC / AUC_WT,
    Basal_norm = Basal / Basal_WT,
    t_MAX_norm = t_MAX / t_MAX_WT
  )

df_norm_summary <- df_norm %>%
  group_by(ID, ID_short, Genotype, Day) %>%
  summarise(
    Max_amplitude_norm = mean(Max_amplitude_norm, na.rm = TRUE),
    AUC_norm = mean(AUC_norm, na.rm = TRUE),
    Basal_norm = mean(Basal_norm, na.rm = TRUE),
    t_MAX_norm = mean(t_MAX_norm, na.rm = TRUE),
    .groups = "drop"
  )

fit <- aov(t_MAX_norm ~ Genotype, data = df_norm_summary)
hist(df_norm$t_MAX_norm)
shapiro.test(resid(fit))
plot(fit)
leveneTest(fit)
t.test(t_MAX_norm ~ Genotype, data = df_norm_summary, var.equal = TRUE)
wilcox.test(df_norm_summary$t_MAX_norm ~ df_norm_summary$Genotype)
