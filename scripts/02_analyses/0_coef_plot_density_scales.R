#####################################################################
#####################################################################
#####################################################################
###
### 0. Testing relationship beween urban density across scales
### Code by: Merin Reji Chacko
### Last edited: 21.01.2026
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list=ls())

explanatory<-read.table("raw_data/explanatory_variables.txt", header = TRUE, sep = ";")

# Pairwise correlations + ggplot heatmap with values and legend forced to 0–1

library(Hmisc)
library(dplyr)
library(tidyr)
library(ggplot2)

density <- explanatory[,c("Urban_50", "Urban_100", "Urban_250", "Urban_500")]

# correlations + p-values
res <- rcorr(as.matrix(density), type = "pearson")
R <- res$r
P <- res$P

# tidy to long form, keep upper triangle only
df <- as.data.frame(as.table(R)) |>
  rename(var1 = Var1, var2 = Var2, r = Freq) |>
  left_join(
    as.data.frame(as.table(P)) |>
      rename(var1 = Var1, var2 = Var2, p = Freq),
    by = c("var1", "var2")
  )

lev <- colnames(density)
df <- df |>
  mutate(
    var1 = factor(var1, levels = lev),
    var2 = factor(var2, levels = lev)
  ) |>
  filter(as.integer(var1) < as.integer(var2))  # upper triangle only

label_fun <- function(x) gsub("_", " ", x)


ggplot(df, aes(x = var2, y = var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r)),
            size = 10, colour = "white") +
  scale_fill_gradient(
    low = "#0077b6",
    high = "#03045e",
    limits = range(df$r, na.rm = TRUE),
    name = "Pearson r"
  )+
  scale_x_discrete(labels = label_fun, position = "bottom") +
  scale_y_discrete(labels = label_fun, position = "right") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, angle = 90, hjust = 0.5, color = "black"),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 11),
    panel.grid = element_blank()
  )

