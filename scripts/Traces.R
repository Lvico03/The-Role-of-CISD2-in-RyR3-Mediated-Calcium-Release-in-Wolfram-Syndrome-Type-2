
library(readr)
library(dplyr)
library(tidyr)
library(purrr) 

files <- list.files(pattern = "\\.csv$", full.names = TRUE)

read_and_transform <- function(file) {
  df <- read_csv(file)
  
  descr_cols <- df %>% select(-starts_with("Mean"))
  mean_cols <- df %>% select(starts_with("Mean"))
  
  bind_cols(descr_cols, mean_cols) %>%
    pivot_longer(
      cols = starts_with("Mean"),
      names_to = "Mean",
      values_to = "Value"
    )
}

long_df <- map_dfr(files, read_and_transform)

write.csv(long_df, "20250603_Hek_RyR3_CISD2_KO_T_Caff.csv", row.names = FALSE)

#Traces 
library(Rmisc)
library(ggsci)
library(tidyverse)
library(splines)

files <- list.files(pattern = "\\.csv$", full.names = TRUE)

dfs <- files %>%
  map(~ read_csv(.x) %>%
        dplyr::rename(Time = `Time (s)`)
  )

df_all <- bind_rows(dfs) %>%
  filter(ID != "20250610_KO_1_4") %>%
  group_by(ID) %>%  
  mutate(Time = as.numeric(Time) - first(as.numeric(Time)) + 100) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = c("RyR3_WT", "RyR3_CISD2_KO")))

p0 <- ggplot(df_all, aes(x = Time, y = Value, group = Genotype, color = Genotype)) +
  geom_smooth(aes(fill = Genotype), method = "lm", alpha = 0.2,
              formula = y ~ splines::bs(x, 20)) +
  scale_color_manual(values = c("RyR3_WT" = "grey25","RyR3_CISD2_KO" = "green4")) +
  scale_fill_manual(values = c("RyR3_WT" = "grey25", "RyR3_CISD2_KO" = "green4")) +
  theme_light() +
  xlab("Time (s)") +
  ylab("FURA-2 (Ex340 / Ex380)") 

print(p0)

library(doBy)
se<-function(x) sd(x) /sqrt(length(x))

#CI use t-distribution rather than z-distribution when deviating from normality and n<30, is what Rmisc CI function does
library(Rmisc)
summary <- summaryBy(Value ~ Time + Genotype, data = df_all,
                     FUN = c(mean = mean, ci = CI))
names(summary)[names(summary) == "Value.FUN1"] <- "Value.mean"   
names(summary)[names(summary) == "Value.FUN2"] <- "CI.mean"     
names(summary)[names(summary) == "Value.FUN3"] <- "CI.lower"     
names(summary)[names(summary) == "Value.FUN4"] <- "CI.upper"

p1 <- ggplot(summary, aes(x = Time, y = Value.mean, color = Genotype, fill = Genotype)) +
  geom_smooth(method = "loess", span = 0.1, se = TRUE, linewidth = 1, alpha = 0.2) +
  scale_color_manual(values = c("RyR3_WT" = "grey25","RyR3_CISD2_KO" = "green4")) +
  scale_fill_manual(values = c("RyR3_WT" = "grey25", "RyR3_CISD2_KO" = "green4")) +
  theme_bw(base_size = 14) +
  xlab("Time (s)") +
  ylab("FURA-2 (Ex340 / Ex380)") 

print(p1)


#Traces derivative 

spar_spline <- 0.5  

# Computation of the derivative 
summary_deriv <- summary %>%
  group_by(Genotype) %>%
  arrange(Time) %>%
  mutate(
    Value.smooth = predict(smooth.spline(Time, Value.mean, spar = spar_spline), Time)$y,
    Deriv_smooth = c(NA, diff(Value.smooth)/diff(Time))
  )

# Traces of the mean derivatives

p_deriv <- ggplot(summary_deriv, aes(x = Time, y = Deriv_smooth, color = Genotype)) +
  geom_line(linewidth = 1.2) +   
  scale_color_manual(values = c("RyR3_WT" = "grey25", "RyR3_CISD2_KO" = "green4")) +
  theme_bw(base_size = 14) +
  xlab("Time (s)") +
  ylab("Mean derivative (ΔValue/ΔTime)")

print(p_deriv)

#Quantification of the slope 

threshold <- 0.01  # considering when the derivative is > 0 -> during the rise of the traces
pendenze <- summary_deriv %>%
  group_by(Genotype) %>%
  filter(Deriv_smooth > threshold) %>%
  summarise(Pendenza_media = mean(Deriv_smooth, na.rm = TRUE))

ggplot(pendenze, aes(x = Genotype, y = Pendenza_media, fill = Genotype)) +
  geom_col(width = 0.4) +
  scale_fill_manual(values = c("RyR3_WT" = "grey25", "RyR3_CISD2_KO" = "green4")) +
  theme_bw(base_size = 16) +
  ylab("Average slope") +
  xlab("Genotype")

