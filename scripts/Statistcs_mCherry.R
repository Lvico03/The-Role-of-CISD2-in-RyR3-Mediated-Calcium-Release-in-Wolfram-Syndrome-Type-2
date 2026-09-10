
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

df_summary <- df_all %>%
  group_by(ID_short, Genotype, Day) %>%
  summarise(
    Max_amplitude = mean(Max_amplitude, na.rm = TRUE),
    AUC = mean(AUC, na.rm = TRUE),
    t_MAX = mean(t_MAX, na.rm = TRUE),
    Basal = mean(Basal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Genotype = factor (Genotype, levels = c("RyR3_WT", "RyR3_Empty", "RyR3_CISD2_CatDead", "RyR3_CISD2_Mutated")))

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

df_norm_summary <- df_norm %>%
  group_by(ID_short, Genotype, Day) %>%
  summarise(
    Max_amplitude_norm = mean(Max_amplitude_norm, na.rm = TRUE),
    AUC_norm = mean(AUC_norm, na.rm = TRUE),
    Basal_norm = mean(Basal_norm, na.rm = TRUE),
    t_MAX_norm = mean(t_MAX_norm, na.rm = TRUE),
    .groups = "drop"
  )

fit <- aov(t_MAX_norm ~ Genotype, data = df_norm_summary)
hist(df_norm_summary$t_MAX_norm)
shapiro.test(resid(fit))
leveneTest(t_MAX_norm ~ Genotype, data = df_norm_summary)

#If data normally distributed 
TukeyHSD(fit)

#If data not normally distributed 
install.packages("FSA")        
library(FSA)
kruskal.test(t_MAX_norm ~ Genotype, data = df_norm_summary)
dunnTest(t_MAX_norm ~ Genotype, data = df_norm_summary, method = "bonferroni")
dunnTest(t_MAX_norm ~ Genotype, data = df_norm_summary, method = "holm")










