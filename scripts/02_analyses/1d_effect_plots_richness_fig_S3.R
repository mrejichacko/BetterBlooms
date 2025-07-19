#####################################################################
#####################################################################
#####################################################################
###
### 1d. Pollinator experiment: Fig S3: 
### Effect plots for richness~ urban_500 
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

## TWO OPTIONS: 
## 1. RUN SCRIPTS 10a-10b FIRST
## 2. LOAD THE PRESAVED ENVIRONMENT as follows:

rm(list=ls())
gc()

load("environments/1a_1b_environment.RData")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# Define pollinator groups with correct sequential labels
pollinator_groups_richness <- list(
  list(label = "(a)", name = "all bees", model = mod.11, bsim = bsim.mod.11, response = "S_allBees"),
  list(label = "(b)", name = "all bees without honeybees", model = mod.12, bsim = bsim.mod.12, response = "S_allBees_noApis"),
  list(label = "(c)", name = "bumblebees", model = mod.14, bsim = bsim.mod.14, response = "S_Bombus"),
  list(label = "(d)", name = "wild social bees", model = mod.15, bsim = bsim.mod.15, response = "S_socialBees"),
  list(label = "(e)", name = "solitary bees", model = mod.16, bsim = bsim.mod.16, response = "S_solitaryBees"),
  list(label = "(f)", name = "wasps", model = mod.17, bsim = bsim.mod.17, response = "S_otherAculeata"),
  list(label = "(g)", name = "hoverflies", model = mod.18, bsim = bsim.mod.18, response = "S_Syrphidae"),
  list(label = "(h)", name = "beetles", model = mod.19, bsim = bsim.mod.19, response = "S_Coleoptera"),
  list(label = "(i)", name = "all pollinators", model = mod.20, bsim = bsim.mod.20, response = "S_allPollinators")
)

# model <- mod.17
# bsim_model <- bsim.mod.17
# response_var <- "S_otherAculeata"

create_prediction_data <- function(model, bsim_model, response_var) {
  newdat <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
  newdat$PlantS.1.z <- quantile(df$PlantS.1.z, probs = c(0.8))  # High plant richness
  newdat$Urban_500.1.z <- (newdat$Urban_500 - mean(df$Urban_500)) / sd(df$Urban_500)
  
  newdat2 <- newdat
  newdat2$PlantS.1.z <- quantile(df$PlantS.1.z, probs = c(0.2))  # Low plant richness
  
  # Design matrix
  Xmat <- model.matrix(~ PlantS.1.z + Urban_500.1.z + PlantS.1.z:Urban_500.1.z, data = newdat)
  
  # Handle different model types
  if ("glmmTMB" %in% class(model)) {
    b <- fixef(model)$cond  # Fixed effects from glmmTMB
  } else if ("glmerMod" %in% class(model)) {
    b <- fixef(model)  # Fixed effects from lme4 (glmer)
  } else {
    b <- coef(model)  # Default for glm and lm
  }
  
  newdat$fit <- exp(Xmat %*% b)  # Point estimates
  
  # Initialize fitmat for confidence intervals
  fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
  
  # Handle different types of Bayesian simulation objects
  if ("sim.merMod" %in% class(bsim_model)) {
    for (i in 1:nsim) {
      b_sim <- as.numeric(bsim_model@fixef[i, ])  # Extract fixed effect samples
      fitmat[, i] <- exp(Xmat %*% b_sim)
    }
  } else {
    for (i in 1:nsim) {
      b_sim <- as.numeric(bsim_model@coef[i, ])  # Standard Bayesian model
      fitmat[, i] <- exp(Xmat %*% b_sim)
    }
  }
  
  newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
  newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)
  
  # Repeat for low plant richness
  Xmat2 <- model.matrix(~ PlantS.1.z + Urban_500.1.z + PlantS.1.z:Urban_500.1.z, data = newdat2)
  newdat2$fit <- exp(Xmat2 %*% b)
  
  fitmat2 <- matrix(ncol = nsim, nrow = nrow(newdat2))
  if ("sim.merMod" %in% class(bsim_model)) {
    for (i in 1:nsim) {
      b_sim2 <- as.numeric(bsim_model@fixef[i, ])
      fitmat2[, i] <- exp(Xmat2 %*% b_sim2)
    }
  } else {
    for (i in 1:nsim) {
      b_sim2 <- as.numeric(bsim_model@coef[i, ])
      fitmat2[, i] <- exp(Xmat2 %*% b_sim2)
    }
  }
  
  newdat2$lwr <- apply(fitmat2, 1, quantile, prob = 0.025)
  newdat2$upr <- apply(fitmat2, 1, quantile, prob = 0.975)
  
  # Label groups
  newdat$group <- "High Plant Richness"
  newdat2$group <- "Low Plant Richness"
  
  # Combine both datasets
  plot_data <- bind_rows(newdat, newdat2)
  plot_data <- subset(plot_data, plot_data$Urban_500 > 0.15 & plot_data$Urban_500 < 0.85)
  
  return(plot_data)
}

#pollinator_group <- pollinator_groups_richness[[6]]

# Function to generate ggplot for richness
create_ggplot <- function(pollinator_group) {
  plot_data <- create_prediction_data(pollinator_group$model, pollinator_group$bsim, pollinator_group$response)
  
  ggplot() +
    geom_point(data = df, aes(x = Urban_500, y = get(pollinator_group$response) / Total_sampling_effort_d), 
               color = "black", alpha = 1, size = 4) +
    geom_ribbon(data = plot_data, aes(x = Urban_500, ymin = lwr, ymax = upr, fill = group), alpha = 0.3) +
    geom_line(data = plot_data, aes(x = Urban_500, y = fit, color = group), size = 1) +
    scale_fill_manual(values = c("High Plant Richness" = "#1B7837", "Low Plant Richness" = "#A6611A")) +
    scale_color_manual(values = c("High Plant Richness" = "#1B7837", "Low Plant Richness" = "#A6611A")) +
    scale_x_continuous(limits = c(0.15, 0.85)) +
    labs(x = "", y = "", title = paste(pollinator_group$label, pollinator_group$name)) +
    theme_classic(base_size = 16) +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0),
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(color = "black", size = 16),  # Make axis numbers black
      axis.ticks = element_line(color = "black")  # Make tick marks black
    )
}

# Generate all plots for richness
all_plots_richness <- lapply(pollinator_groups_richness, create_ggplot)

# Create a lookup for plot labels
label_lookup <- c(
  "Anthophila"        = "bee",
  "Anthophila_noApis"    = "bee no honey bee",
  "Bombus"            = "bumblebee",
  "social"            = "wild social bee",
  "solitary"          = "solitary bee",
  "other_Aculeata"    = "wasp",
  "Syrphidae"         = "hoverfly",
  "Coleoptera"        = "beetle",
  "All_pollinators"   = "all pollinator"
)

# Define custom colors and plant order
latin_names <- c(
  carrot   = "D. carota",
  radish   = "R. sativus",
  sainfoin = "O. viciifolia",
  comfrey  = "S. officinale",
  all      = "All plants"
)

plant_order <- c("D. carota", "R. sativus", "O. viciifolia", "S. officinale")

custom_colors <- c(
  "D. carota" = "#FCAA67",
  "R. sativus" = "#D33665",
  "O. viciifolia" = "#E08DAC",
  "S. officinale" = "#8D80AD"
)

div <- read.table("cleaned_data/pollinators_sum_stat.txt", header = T, sep = ";")

# Filter out the one group we want to skip
div_filtered <- div %>%
  filter(pollinator_group != "Apis mellifera")

# Create an empty list to store the patchwork plots
plot_list <- list()

# Loop over rows
for (i in 1:nrow(div_filtered)) {
  row <- div_filtered[i, ]
  poll_group <- row$pollinator_group
  label <- label_lookup[poll_group]
  
  # Prepare data
  row_long <- row %>%
    pivot_longer(cols = -pollinator_group, names_to = "metric", values_to = "value") %>%
    separate(metric, into = c("type", "plant"), sep = "_") %>%
    mutate(
      type = ifelse(type == "A", "Abundance", "Richness"),
      plant_latin = recode(plant, !!!latin_names),
      plant_latin = factor(plant_latin, levels = plant_order)
    )
  
  # Abundance data
  abund_data <- row_long %>% filter(type == "Abundance", plant != "all")
  
  # Richness data
  rich_data <- row_long %>% filter(type == "Richness")
  rich_plot_data <- rich_data %>% filter(plant != "all")
  total_richness <- filter(rich_data, plant == "all")$value
  
  # Abundance plot
  abund_plot <- ggplot(abund_data, aes(x = "Total", y = value, fill = plant_latin)) +
    geom_bar(stat = "identity") +
    labs(title = "", y = "Total abundance", x = "") +
    scale_fill_manual(values = custom_colors) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 16, color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.title = element_text(size = 18, color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
      
    )
  
  # Richness plot
  richness_label <- paste(label, "species richness")
  # Richness plot
  rich_plot <- ggplot(rich_plot_data, aes(x = plant_latin, y = value, fill = plant_latin)) +
    geom_bar(stat = "identity", alpha = 1) +
    geom_hline(yintercept = total_richness, color = "black", linetype = "solid") +
    annotate("text", x = 3, y = total_richness + 3,
             label = "Total richness", hjust = 0.5, size = 4) +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = "",
      y = "Richness",
      x = ""
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 16, color = "black", 
                                 angle = 90, hjust = 1,
                                 face = "italic"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18, color = "black")
    )
  
  if (!is.na(total_richness) && length(total_richness) > 0 && total_richness > 0) {
    rich_plot <- rich_plot +
      geom_hline(yintercept = total_richness, color = "black", linetype = "solid")
  }
  
  # Combine with model plot (assumed to be all_plots[[i]]) — make sure you align i with the loop index!
  model_plot <- all_plots_richness[[i]]
  model_plot <- model_plot +
    labs(y = "Pollinator richness\n(spp. / 9h sampling)", x = "Densification")
  
  # Combine model plot and abundance plot only
  final_plot <- (model_plot | rich_plot) + plot_layout(widths = c(4, 1))
  
  # Store it
  plot_list[[poll_group]] <- final_plot
  
}

# Show one:
library(patchwork)

wrap_plots(plot_list, ncol = 3)

# save landscape 24 x 20
