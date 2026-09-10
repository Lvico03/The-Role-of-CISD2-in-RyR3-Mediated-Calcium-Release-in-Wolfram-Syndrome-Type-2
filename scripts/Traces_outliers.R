
library(Rmisc)
library(tidyverse)
library(ggsci)

files <- list.files(pattern = "\\.csv$", full.names = TRUE)

dfs <- files %>%
  map(~ read_csv(.x) %>%
        rename(Time = `Time (s)`) %>%
        mutate(Time = as.numeric(Time) - as.numeric(Time[1]) + 100)
  )

df_all <- bind_rows(dfs) %>%
  filter(ID != "20250617_CISD2_KO_1_1") %>% #outlier well removal
  mutate(Genotype = factor(Genotype, levels = c("RyR3_WT", "RyR3_CISD2_KO")))
  
  mean_WT <- df_all %>%
    filter(Genotype == "RyR3_WT") %>%
    summarise(mean_WT = mean(Value, na.rm = TRUE)) %>%
    pull(mean_WT)
  
  df_all <- df_all %>%
    mutate(Value_norm = Value / mean_WT)
  
names(read_csv(files[1])) 

# Outlaiers
outlier_summary <- df_all %>%
  group_by(ID) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val = sd(Value, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(mean_val))

print(outlier_summary)

p_box <- ggplot(df_all, aes(x = ID, y = Value)) +
  geom_boxplot() +
  coord_flip() +
  theme_light() +
  ggtitle("ID values distribution")

print(p_box)

Q1 <- quantile(df_all$Value, 0.25, na.rm = TRUE)
Q3 <- quantile(df_all$Value, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

outliers <- df_all %>%
  filter(Value < (Q1 - 1.5 * IQR) | Value > (Q3 + 1.5 * IQR))

p_traces <- ggplot(df_all, aes(x = Time, y = Value, group = ID, color = ID)) +
  geom_line(alpha = 0.3) +
  geom_point(data = outliers, aes(x = Time, y = Value), color = "red", size = 2) +
  theme_light() +
  ggtitle("Traces with outliers")

print(p_traces)

id_counts <- df_all %>% count(ID)
print(id_counts)

#oulier removal
df_filtered <- df_all %>%
  filter(ID != "20250617_CISD2_KO_1_1")

p0 <- ggplot(df_all, aes(x = Time, y = Value_norm, group = Genotype, color = Genotype)) +
  geom_smooth(aes(fill = Genotype), method = "lm", alpha = 0.2,
              formula = y ~ splines::bs(x, 20)) +
  theme_light() +
  xlab("Time (s)") +
  ylab("Normalized FURA-2 (Ex340 / Ex380)") +
  xlim(c(80, 560)) +
  ggtitle("Smoothed traces (normalized to WT, outlier removed)")

print(p0)

# Mean ± CI (normalized) + Time binning
df_all <- df_all %>%
  mutate(Time_bin = floor(Time / 10) * 10)

summary_df <- summarySE(df_all, measurevar = "Value_norm", groupvars = c("Time_bin", "Genotype"))

p1 <- ggplot(summary_df, aes(x = Time_bin, y = Value_norm, group = Genotype)) +
  geom_path(aes(colour = Genotype), linewidth = 1) +
  geom_ribbon(aes(ymin = Value_norm - ci, ymax = Value_norm + ci, fill = Genotype), alpha = 0.2, linetype = "blank") +
  scale_color_manual(values = c("RyR3_WT" = "grey3", "RyR3_CISD2_KO" = "green4")) +
  scale_fill_manual(values = c("RyR3_WT" = "grey3", "RyR3_CISD2_KO" = "green4")) +
  theme_bw(base_size = 14) +
  xlab("Time (s, binned every 10s)") +
  ylab("Normalized FURA-2 (Ex340 / Ex380)") +
  xlim(c(80, 560)) +
  ggtitle("Binned means ± CI (normalized to WT, outlier removed)")

print(p1)
