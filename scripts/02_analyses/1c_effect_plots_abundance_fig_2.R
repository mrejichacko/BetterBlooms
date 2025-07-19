#####################################################################
#####################################################################
#####################################################################
###
### 1c. Pollinator experiment: figure 2
### Effect plots for abundance ~ urban_500 
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

## TWO OPTIONS: 
## 1. RUN SCRIPTS 1a-1b FIRST, OR
## 2. LOAD THE PRESAVED ENVIRONMENT as follows:

rm(list=ls())
gc()

load("environments/1a_1b_environment.RData")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(patchwork)
library(tidyverse)

# Define pollinator groups with correct sequential labels
pollinator_groups <- list(
  list(label = "(a)", name = "all bees", model = mod.1, bsim = bsim.mod.1, response = "A_allBees"),
  list(label = "(b)", name = "honeybees", model = mod.3, bsim = bsim.mod.3, response = "A_Apis"),
  list(label = "(c)", name = "bumblebees", model = mod.4, bsim = bsim.mod.4, response = "A_Bombus"),
  list(label = "(d)", name = "wild social bees", model = mod.5, bsim = bsim.mod.5, response = "A_socialBees"),
  list(label = "(e)", name = "solitary bees", model = mod.6, bsim = bsim.mod.6, response = "A_solitaryBees"),
  list(label = "(f)", name = "wasps", model = mod.7, bsim = bsim.mod.7, response = "A_otherAculeata"),
  list(label = "(g)", name = "hoverflies", model = mod.8, bsim = bsim.mod.8, response = "A_Syrphidae"),
  list(label = "(h)", name = "beetles", model = mod.9, bsim = bsim.mod.9, response = "A_Coleoptera"),
  list(label = "(i)", name = "all pollinators", model = mod.10, bsim = bsim.mod.10, response = "A_allPollinators")
)



# Function to create prediction dataset
create_prediction_data <- function(model, bsim_model, response_var) {
  newdat <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
  newdat$PlantS.1.z <- quantile(df$PlantS.1.z, probs = c(0.8))  # High plant richness
  newdat$Urban_500.1.z <- (newdat$Urban_500 - mean(df$Urban_500)) / sd(df$Urban_500)
  
  newdat2 <- newdat
  newdat2$PlantS.1.z <- quantile(df$PlantS.1.z, probs = c(0.2))  # Low plant richness
  
  # Design matrix
  Xmat <- model.matrix(~ PlantS.1.z + Urban_500.1.z + PlantS.1.z:Urban_500.1.z, data=newdat)
  b <- fixef(model)
  newdat$fit <- exp(Xmat %*% b)
  
  # Simulate confidence intervals using precomputed bsim model
  fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
  for (i in 1:nsim) fitmat[, i] <- exp(Xmat %*% bsim_model@fixef[i,])
  newdat$lwr <- apply(fitmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(fitmat, 1, quantile, prob=0.975)
  
  # Repeat for low plant richness
  Xmat2 <- model.matrix(~ PlantS.1.z + Urban_500.1.z + PlantS.1.z:Urban_500.1.z, data=newdat2)
  newdat2$fit <- exp(Xmat2 %*% b)
  
  fitmat2 <- matrix(ncol = nsim, nrow = nrow(newdat2))
  for (i in 1:nsim) fitmat2[, i] <- exp(Xmat2 %*% bsim_model@fixef[i,])
  newdat2$lwr <- apply(fitmat2, 1, quantile, prob=0.025)
  newdat2$upr <- apply(fitmat2, 1, quantile, prob=0.975)
  
  # Label the groups
  newdat$group <- "High Plant Richness"
  newdat2$group <- "Low Plant Richness"
  
  # Combine both datasets
  plot_data <- bind_rows(newdat, newdat2)
  plot_data <- subset(plot_data, plot_data$Urban_500 > 0.15 & plot_data$Urban_500 < 0.85)
  
  return(plot_data)
}

# Modify the plot title inside `create_ggplot()` function:
create_ggplot <- function(pollinator_group) {
  plot_data <- create_prediction_data(pollinator_group$model, pollinator_group$bsim, pollinator_group$response)
  
  ggplot() +
    geom_point(data = df, aes(x = Urban_500, y = get(pollinator_group$response) / Total_sampling_effort_d), 
               color = "black", alpha = 1, size = 4) +
    geom_ribbon(data = plot_data, aes(x = Urban_500, ymin = lwr, ymax = upr, fill = group), alpha = 0.3) +
    geom_line(data = plot_data, aes(x = Urban_500, y = fit, color = group), linewidth = 1) +
    scale_fill_manual(values = c("High Plant Richness" = "#1B7837", "Low Plant Richness" = "#A6611A")) +
    scale_color_manual(values = c("High Plant Richness" = "#1B7837", "Low Plant Richness" = "#A6611A")) +
    scale_x_continuous(limits = c(0.15, 0.85)) +
    labs(x = "", y = "", title = paste(pollinator_group$label, pollinator_group$name)) +
    theme_classic(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0),
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(color = "black", size = 16),  # Make axis numbers black
      axis.ticks = element_line(color = "black")  # Make tick marks black
    )
  
}

# Generate all plots
all_plots <- lapply(pollinator_groups, create_ggplot)

## ADD THE ABUNDANCE DATA

# Create a lookup for plot labels
label_lookup <- c(
  "Anthophila"        = "bee",
  "Apis mellifera"    = "honey bee",
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
  carrot   = "Daucus\ncarota",
  radish   = "Raphanus\nsativus",
  sainfoin = "Onobrychis\nviciifolia",
  comfrey  = "Symphytum\nofficinale",
  all      = "All plants"
)

plant_order <- c("Daucus\ncarota", "Raphanus\nsativus", "Onobrychis\nviciifolia", "Symphytum\nofficinale")

custom_colors <- c(
  "Daucus\ncarota" = "#FCAA67",
  "Raphanus\nsativus" = "#D33665",
  "Onobrychis\nviciifolia" = "#E08DAC",
  "Symphytum\nofficinale" = "#8D80AD"
)

div <- read.table("cleaned_data/pollinators_sum_stat.txt", header = T, sep = ";")
# Filter out the one group we want to skip
div_filtered <- div %>%
  filter(pollinator_group != "Anthophila_noApis")

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
  
  # Optional total richness line (omit for now unless you want to include it for certain groups)
  if (!is.na(total_richness) && length(total_richness) > 0 && total_richness > 0) {
    rich_plot <- rich_plot +
      geom_hline(yintercept = total_richness, color = "black", linetype = "solid")
  }
  
  # Combine with model plot (assumed to be all_plots[[i]]) — make sure you align i with the loop index!
  model_plot <- all_plots[[i]]
  model_plot <- model_plot +
    labs(y = "Flower visitor frequency\n(indiv. / 9h sampling)", x = "Proportion of impervious surface\n(buffer of 500-m around garden)")
  
  # Combine model plot and abundance plot only
  final_plot <- (model_plot | abund_plot) + plot_layout(widths = c(4, 1))
  
  # Store it
  plot_list[[poll_group]] <- final_plot
  
  
  # Store combined layout
  plot_list[[poll_group]] <- final_plot
  
}


# Show one:
library(patchwork)

wrap_plots(plot_list, ncol = 3)

# save landscape 25 x 20
