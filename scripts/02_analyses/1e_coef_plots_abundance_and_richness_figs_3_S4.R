#####################################################################
###
### 1e. Pollinator experiment: abundance, plotting coefs
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list=ls())

library(tidyverse)
library(ggplot2)

dd <- read.table("results/beta_coefficients_all_plants_for_paper.txt",
                 sep = ";", header = T)
dd$Model
dd$Model_name <- c("All bees", "Bees without honeybees", "Honeybees", "Bumblebees", "Social wild bees", "Solitary bees","Wasps", "Hoverflies", "Beetles", "All pollinators",
                   "All bees", "Bees without honeybees", "Bumblebees", "Social wild bees", "Solitary bees","Wasps", "Hoverflies", "Beetles", "All pollinators")
dd$N <- 23

dd$Type <- ""
dd$Type[1:10] <- "Abundance"
dd$Type[11:19] <- "Richness"
names(dd) <- c("Model",
               "Mean_Intercept",
               "2.5_CI_Intercept",
               "97.5_CI_Intercept",
               "Mean_Local",
               "2.5_CI_Local",
               "97.5_CI_Local",
               "Mean_Landscape",
               "2.5_CI_Landscape",
               "97.5_CI_Landscape",
               "Mean_Interaction",
               "2.5_CI_Interaction",
               "97.5_CI_Interaction",
               "Model_name",
               "N",
               "Type")


# Filter abundance only, pivot longer
dd_abund <- dd %>%
  filter(Type == "Abundance") %>%
  select(Model_name, starts_with("Mean_"), starts_with("2.5_CI_"), starts_with("97.5_CI_")) %>%
  pivot_longer(
    cols = -Model_name,
    names_to = c(".value", "Predictor"),
    names_pattern = "(Mean|2.5_CI|97.5_CI)_(.*)"
  ) %>%
  mutate(
    Predictor = factor(Predictor, levels = c("Local", "Landscape", "Interaction")),
    Sig = ifelse(`2.5_CI` > 0 | `97.5_CI` < 0, "Credible", "Uncertain")
  )

dd_abund <- subset(dd_abund, !is.na(dd_abund$Predictor))
dd_abund$Predictor = factor(dd_abund$Predictor, levels = c("Interaction", "Landscape", "Local"))

custom_labels <- c(
  "(a) all bees",
  "bees without honeybees",
  "(b) honeybees",
  "(c) bumblebees",
  "(d) social wild bees",
  "(e) solitary bees",
  "(f) wasps",
  "(g) hoverflies",
  "(h) beetles",
  "(i) all pollinators"
)

# Recode Model_name for abundance models only
dd_abund <- dd_abund %>%
  mutate(Model_name = factor(Model_name,
                             levels = c("All bees",
                                        "Bees without honeybees",
                                        "Honeybees",
                                        "Bumblebees",
                                        "Social wild bees",
                                        "Solitary bees",
                                        "Wasps",
                                        "Hoverflies",
                                        "Beetles",
                                        "All pollinators"),
                             labels = custom_labels))

dd_abund <- subset(dd_abund, Model_name != "bees without honeybees")

ggplot(dd_abund, aes(x = Mean, y = Predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(aes(color = Sig), size = 4) +
  geom_errorbarh(aes(xmin = `2.5_CI`, xmax = `97.5_CI`, color = Sig), height = 0, linewidth=0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  facet_wrap(~Model_name, ncol = 3, scales = "free_x") +
  labs(
    x = "Mean parameter estimate and 95% credible interval",
    y = "",
    title = "",
    color = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",  # Keep legend for shape
    strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black"),
    #axis.ticks.y= element_blank()
  )

#for now: saved 10 x 10 inches

### richness

dd_rich <- dd %>%
  filter(Type == "Richness") %>%
  select(Model_name, starts_with("Mean_"), starts_with("2.5_CI_"), starts_with("97.5_CI_")) %>%
  pivot_longer(
    cols = -Model_name,
    names_to = c(".value", "Predictor"),
    names_pattern = "(Mean|2.5_CI|97.5_CI)_(.*)"
  ) %>%
  mutate(
    Predictor = factor(Predictor, levels = c("Interaction", "Landscape","Local")),
    Sig = ifelse(`2.5_CI` > 0 | `97.5_CI` < 0, "Credible", "Uncertain")
  ) %>%
  filter(!is.na(Predictor))

custom_labels_richness <- c(
  "(a) all bees",
  "(b) bees without honeybees",
  "(c) bumblebees",
  "(d) social wild bees",
  "(e) solitary bees",
  "(f) wasps",
  "(g) hoverflies",
  "(h) beetles",
  "(i) all pollinators"
)

dd_rich <- dd_rich %>%
  mutate(Model_name = factor(Model_name,
                             levels = c("All bees",
                                        "Bees without honeybees",
                                        "Bumblebees",
                                        "Social wild bees",
                                        "Solitary bees",
                                        "Wasps",
                                        "Hoverflies",
                                        "Beetles",
                                        "All pollinators"),
                             labels = custom_labels_richness))

ggplot(dd_rich, aes(x = Mean, y = Predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(aes(color = Sig), size = 4) +
  geom_errorbarh(aes(xmin = `2.5_CI`, xmax = `97.5_CI`, color = Sig), height = 0, linewidth = 0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  facet_wrap(~Model_name, ncol = 3, scales = "free_x") +
  labs(
    x = "Mean parameter estimate and 95% credible interval",
    y = "",
    title = "",
    color = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black")
  )
# save 10 x 10 for now, the rest of the information was added in inkscape