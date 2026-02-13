#####################################################################
#####################################################################
#####################################################################
###
### 1a.2. Modelling abundance by phytometer
### Code by: Merin Reji Chacko
### Last edited: 14.01.2026
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

rm(list=ls())

floristic<-read.table("raw_data/plant_floristic_data.txt", header = TRUE, sep = ";")
explanatory<-read.table("raw_data/explanatory_variables.txt", header = TRUE, sep = ";")
abundance<-read.table("cleaned_data/pollinator_abundance_aggregated_garden_phytometer.txt", header = TRUE, sep = ";")
richness<-read.table("cleaned_data/pollinator_SR_aggregated_garden_phytometer.txt", header = TRUE, sep = ";")
effort<-read.table("cleaned_data/sampling_effort_gardens_plants_aggregated.txt", header = TRUE, sep = ";")

str(floristic)
str(explanatory)
str(abundance)
str(richness)
str(effort)

dim(explanatory)

#Combine the datasets
df <- Reduce(function(x, y) merge(x, y, by="Id", sort=TRUE), 
             list(floristic[,c(1,5)],
                  explanatory[,c(1,4,5,13:16)],
                  abundance,
                  richness,
                  effort))

summary(df)
dim(df)
str(df)
names(df)

#Rename the plant variable
colnames(df)[2]<-"PlantS"
colnames(df)[2]

#Kick out garden Nr. 39, the sampling effort is too low
df <- droplevels(df[which(df$Id != 39),])
dim(df)

rm(abundance, effort, explanatory, floristic, richness)

###########################
### Scale the variables ###
###########################

df$Id.fac <- as.factor(df$Id)

##Landscape, all scales:

poly.Urban_500 <- poly(df$Urban_500,2)
df$Urban_500.1.z <- scale(poly.Urban_500[,1])
df$Urban_500.2.z <- scale(poly.Urban_500[,2])

poly.Urban_250 <- poly(df$Urban_250,2)
df$Urban_250.1.z <- scale(poly.Urban_250[,1])
df$Urban_250.2.z <- scale(poly.Urban_250[,2])

poly.Urban_100 <- poly(df$Urban_100,2)
df$Urban_100.1.z <- scale(poly.Urban_100[,1])
df$Urban_100.2.z <- scale(poly.Urban_100[,2])

poly.Urban_50 <- poly(df$Urban_50,2)
df$Urban_50.1.z <- scale(poly.Urban_50[,1])
df$Urban_50.2.z <- scale(poly.Urban_50[,2])

##Garden:

#Plants species richness May-August
poly.PlantS <- poly(df$PlantS, 2)
df$PlantS.1.z <- scale(poly.PlantS[,1])
df$PlantS.2.z <- scale(poly.PlantS[,2])

##Sampling:

#Check the sampling variables and transform them in daily rates:
df$Total_sampling_effort_d <- df$total_sampling_effort_min/60/9 #first minutes to h, then  per fieldwork day = 9h (one fieldwork day lasted from 9:00 to 18:00)
hist(df$Total_sampling_effort_d) 

df$n_flowers_carrot.z<-scale(df$n_flowers_carrot)
df$n_flowers_radish.z<-scale(df$n_flowers_radish)
df$n_flowers_sainfoin.z<-scale(df$n_flowers_sainfoin)
df$n_flowers_comfrey.z<-scale(df$n_flowers_comfrey)

df$n_flowers.z <-scale((df$n_flowers_carrot+df$n_flowers_radish+df$n_flowers_sainfoin+df$n_flowers_comfrey))

## -----------------------------
## Per-phytometer abundance GLMMs
## -----------------------------

library(lme4)
library(blmeco)
library(plyr)
library(AICcmodavg)
library(arm)

nsim <- 10000

phyto_levels <- c("Carrot","Radish","Sainfoin","Comfrey")

poll_groups <- c("allBees","allBees_noApis","Apis","Bombus","socialBees",
                 "solitaryBees","otherAculeata","Syrphidae","Coleoptera","allPollinators")

resp_map <- setNames(
  c("A_allBees","A_allBees_noApis","A_Apis","A_Bombus","A_socialBees",
    "A_solitaryBees","A_otherAculeata","A_Syrphidae","A_Coleoptera","A_allPollinators"),
  poll_groups
)

effort_map <- c(
  Carrot   = "sampling_effort_min_carrot",
  Radish   = "sampling_effort_min_radish",
  Sainfoin = "sampling_effort_min_sainfoin",
  Comfrey  = "sampling_effort_min_comfrey"
)

## build effort per day per phyto
for (p in phyto_levels) {
  evar <- effort_map[[p]]
  newe <- paste0("effort_d_", p)
  df[[newe]] <- df[[evar]]/60/9
  df[[newe]] <- ifelse(df[[newe]] > 0, df[[newe]], NA_real_)
}

fit_one_phyto <- function(p) {
  
  effort_var <- paste0("effort_d_", p)
  
  mods <- list()
  
  for (g in poll_groups) {
    resp <- paste0(resp_map[[g]], "_", p)
    
    if (!resp %in% names(df)) next
    
    dsub <- df[!is.na(df[[effort_var]]), , drop = FALSE]
    
    f <- as.formula(
      paste0(resp, " ~ PlantS.1.z * Urban_500.1.z + (1|Id.fac) + offset(log(", effort_var, "))")
    )
    
    mods[[g]] <- glmer(f, data = dsub, family = poisson())
  }
  
  ## model fit table
  fit_tab <- ldply(mods, function(m) {
    data.frame(
      AICc = round(AICc(m), 2),
      DispPar = dispersion_glmer(m),
      n = nobs(m)
    )
  })
  colnames(fit_tab)[1] <- "Pollinator_group"
  fit_tab$Phytometer <- p
  
  ## coefficient intervals via arm::sim
  coef_tab <- ldply(mods, function(m) {
    simres <- sim(m, n.sim = nsim)@fixef
    q <- apply(simres, 2, quantile, prob = c(0.025, 0.5, 0.975))
    out <- as.data.frame(t(q))
    out$fixeff <- rownames(out)
    out
  })
  colnames(coef_tab)[1] <- "Pollinator_group"
  colnames(coef_tab)[2:4] <- c("lower","mean","upper")
  coef_tab$Phytometer <- p
  coef_tab$Response <- "Abundance"
  coef_tab$Scale <- 500
  
  list(models = mods, fit = fit_tab, coefs = coef_tab)
}

res_by_phyto <- lapply(phyto_levels, fit_one_phyto)
names(res_by_phyto) <- phyto_levels

fit_all  <- do.call(rbind, lapply(res_by_phyto, `[[`, "fit"))
coefs_all <- do.call(rbind, lapply(res_by_phyto, `[[`, "coefs"))

# write.table(fit_all,
#             file = "results/abundance_models_model_fit_by_phytometer_500m.txt",
#             sep = ";", row.names = FALSE, quote = TRUE)
# 
# write.table(coefs_all,
#             file = "results/abundance_models_coefficients_by_phytometer_500m.txt",
#             sep = ";", row.names = FALSE, quote = TRUE)

library(ggplot2)
library(dplyr)
library(stringr)

coefs <- coefs_all
coefs_pooled <- read.table("results/abundance_models_coefficents_all_scales.txt", header = TRUE, sep = ";")

# Build pooled rows to append to per-phyto table
pooled_as_level <- coefs_pooled %>%
  filter(Scale == 500) %>%
  transmute(
    Pollinator_group,
    fixeff,
    lower,
    mean,
    upper,
    Phytometer = "Pooled"
  )%>%
  mutate(
    Pollinator_group = recode(Pollinator_group,
                              "allBees" = "All bees",
                              "allBees_noApis" = "All bees (no honeybees)",
                              "Apis" = "Honeybees",
                              "Bombus" = "Bumblebees",
                              "socialBees" = "Social bees",
                              "solitaryBees" = "Solitary bees",
                              "otherAculeata" = "Wasps",
                              "Syrphidae" = "Hoverflies",
                              "Coleoptera" = "Beetles",
                              "allPollinators" = "All pollinators"
    ),
    fixeff = factor(fixeff, levels = levels(coefs$fixeff))
  )


# optional: nicer labels
coefs <- coefs %>%
  mutate(
    Pollinator_group = recode(Pollinator_group,
                              "allBees" = "All bees",
                              "allBees_noApis" = "All bees (no honeybees)",
                              "Apis" = "Honeybees",
                              "Bombus" = "Bumblebees",
                              "socialBees" = "Social bees",
                              "solitaryBees" = "Solitary bees",
                              "otherAculeata" = "Wasps",
                              "Syrphidae" = "Hoverflies",
                              "Coleoptera" = "Beetles",
                              "allPollinators" = "All pollinators"
    ),
    Phytometer = factor(Phytometer, levels = c("Carrot","Radish","Sainfoin","Comfrey")),
    fixeff = factor(fixeff, levels = c("(Intercept)","PlantS.1.z","Urban_500.1.z","PlantS.1.z:Urban_500.1.z"))
  )

plot_df <- bind_rows(
  coefs,                # your per-phytometer table
  pooled_as_level
) %>%
  mutate(
    Phytometer = factor(Phytometer, levels = c("Carrot","Radish","Sainfoin","Comfrey","Pooled"))
  )


plot_df <- plot_df %>%
  mutate(
    fixeff = factor(
      fixeff,
      levels = c(
        "(Intercept)",
        "PlantS.1.z",
        "Urban_500.1.z",
        "PlantS.1.z:Urban_500.1.z"
      ),
      labels = c(
        "Intercept",
        "Floral richness",
        "Densification",
        "Floral richness × Densification"
      )
    )
  )


plot_df <- plot_df %>%
  mutate(
    sig_class = case_when(
      lower > 0  ~ "positive",
      upper < 0  ~ "negative",
      TRUE       ~ "nonsignificant"
    ),
    sig_class = factor(
      sig_class,
      levels = c("positive", "negative", "nonsignificant")
    )
  )

group_order <- c(
  "All bees",
  "All bees (no honeybees)",
  "Honeybees",
  "Bumblebees",
  "Social bees",
  "Solitary bees",
  "Wasps",
  "Hoverflies",
  "Beetles",
  "All pollinators"
)

plot_df <- plot_df %>%
  mutate(
    Pollinator_group = factor(Pollinator_group, levels = group_order)
  )

plot_df <- subset(plot_df, Pollinator_group != "All bees (no honeybees)")

labeller_groups <- c(
  "All bees"        = "(a) all bees",
  "Honeybees"      = "(b) honeybees",
  "Bumblebees"     = "(c) bumblebees",
  "Social bees"    = "(d) social bees",
  "Solitary bees"  = "(e) solitary bees",
  "Wasps"          = "(f) wasps",
  "Hoverflies"     = "(g) hoverflies",
  "Beetles"        = "(h) beetles",
  "All pollinators"= "(i) all pollinators"
  
)

phyto_labels <- c(
  "Carrot"   = "italic('D. carota')",
  "Radish"   = "italic('R. sativus')",
  "Sainfoin" = "italic('O. viciifolia')",
  "Comfrey"  = "italic('S. officinale')",
  "Pooled"   = "Pooled"
)



ggplot(plot_df %>% filter(fixeff == "Densification"), aes(x = Phytometer, y = mean, color = sig_class)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  coord_flip()+
  geom_errorbar(linewidth = .5, aes(ymin = lower, ymax = upper), width = 0) +
  facet_wrap(~Pollinator_group, scales = "free_x", nrow=3,
             labeller = labeller(Pollinator_group = labeller_groups)) +
  scale_x_discrete(
    limits = rev(levels(plot_df$Phytometer)),
    labels = function(x) parse(text = phyto_labels[x])
  )+
  scale_colour_manual(
    values = c(
      "positive" = "#669bbc",   # blue
      "negative" = "#c1121f",   # red
      "nonsignificant" = "black"
    ),
    labels = c(
      "positive" = "Positive (95% CI > 0)",
      "negative" = "Negative (95% CI < 0)",
      "nonsignificant" = "Not significant"
    ),
    name = "Effect"
  ) +
  labs(
    x = NULL,
    y = "Mean parameter estimate and 95% credible interval"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    axis.text = element_text(color = "black"),
    legend.position = "bottom")

# portrait 9 x 8 inches    
