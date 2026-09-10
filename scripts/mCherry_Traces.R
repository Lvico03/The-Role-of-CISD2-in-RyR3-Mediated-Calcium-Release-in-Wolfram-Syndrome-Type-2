
#mCherry Traces 
library(Rmisc)
library(ggsci)
library(tidyverse)

files <- list.files(pattern = "\\.csv$", full.names = TRUE)

dfs <- files %>%
  map(~ read_csv(.x) %>%
        rename(Time = `Time (s)`)
  )

df_all <- bind_rows(dfs) %>%
  group_by(ID) %>%  
  mutate(Time = as.numeric(Time) - first(as.numeric(Time)) + 100) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = c("RyR3_WT", "RyR3_Empty", "RyR3_CISD2_CatDead", "RyR3_CISD2_Mutated")))

p0 <- ggplot(df_all, aes(x = Time, y = Value, group = Genotype, color = Genotype)) +
  geom_smooth(aes(fill = Genotype), method = "lm", alpha = 0.2,
              formula = y ~ splines::bs(x, 20)) +
  scale_color_manual(values = c("RyR3_WT" = "grey25","RyR3_Empty" = "green4", "RyR3_CISD2_CatDead" = "blue", "RyR3_CISD2_Mutated" = "purple")) +
  scale_fill_manual(values = c("RyR3_WT" = "grey25", "RyR3_Empty" = "green4", "RyR3_CISD2_CatDead" = "blue", "RyR3_CISD2_Mutated" = "purple")) +
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
  scale_color_manual(values = c("RyR3_WT" = "grey25","RyR3_Empty" = "green4", "RyR3_CISD2_CatDead" = "blue", "RyR3_CISD2_Mutated" = "purple")) +
  scale_fill_manual(values = c("RyR3_WT" = "grey25", "RyR3_Empty" = "green4", "RyR3_CISD2_CatDead" = "blue", "RyR3_CISD2_Mutated" = "purple")) +
  theme_bw(base_size = 14) +
  xlab("Time (s)") +
  ylab("FURA-2 (Ex340 / Ex380)") 

print(p1)

