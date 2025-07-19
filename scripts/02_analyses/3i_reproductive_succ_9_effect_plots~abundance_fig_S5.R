#####################################################################
#####################################################################
#####################################################################
###
### Pollinator experiment: Modelling plant reproductive success: Effect plots for abundance 
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
#####################################################################
#####################################################################
#####################################################################
###
### Prep

#####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
### RUN 3a-h!!! It's not enough to load the 3a environment, you need all the models!! ####
#####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###

library(ggplot2)
#########################
# (a) Carrot: seed set  #
#########################

#Show the effect for the hoverflies
newdat_mod.1_Carrot <- data.frame(A_Syrphidae_Carrot.dayly = seq(0, 18, by = 1))
newdat_mod.1_Carrot$A_Syrphidae_Carrot.dayly.z <- (newdat_mod.1_Carrot$A_Syrphidae_Carrot.dayly-mean(df7a$A_Syrphidae_Carrot.dayly))/sd(df7a$A_Syrphidae_Carrot.dayly)
newdat_mod.1_Carrot$A_Apis_Carrot.dayly.z <- 0 
newdat_mod.1_Carrot$A_socialBees_Carrot.dayly.z <- 0 
newdat_mod.1_Carrot$A_solitaryBees_Carrot.dayly.z <- 0 
newdat_mod.1_Carrot$A_otherAculeata_Carrot.dayly.z <- 0
newdat_mod.1_Carrot$A_Coleoptera_Carrot.dayly.z <-0

Xmat_mod.1_Carrot <- model.matrix(~ A_Apis_Carrot.dayly.z + 
                                    A_socialBees_Carrot.dayly.z + 
                                    A_solitaryBees_Carrot.dayly.z + 
                                    A_otherAculeata_Carrot.dayly.z + 
                                    A_Syrphidae_Carrot.dayly.z + A_Coleoptera_Carrot.dayly.z, data=newdat_mod.1_Carrot)
newdat_mod.1_Carrot$fit <- exp(Xmat_mod.1_Carrot %*% fixef(mod.1_Carrot))# exp = inverse link function for poisson models. 
fitmat_mod.1_Carrot <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Carrot))
for(i in 1:nsim) fitmat_mod.1_Carrot[,i] <- exp(Xmat_mod.1_Carrot %*% bsim.mod.1_Carrot@fixef[i,])
newdat_mod.1_Carrot$lwr <- apply(fitmat_mod.1_Carrot, 1, quantile, prob=0.025)
newdat_mod.1_Carrot$upr <- apply(fitmat_mod.1_Carrot, 1, quantile, prob=0.975)

a <- ggplot(newdat_mod.1_Carrot, aes(x = A_Syrphidae_Carrot.dayly, y = fit)) +
  geom_point(data = df7a, aes(x = A_Syrphidae_Carrot.dayly, y = n_seeds), 
             color = "black", 
             shape = 17, 
             fill = "black", 
             size = 4,
             stroke = 0, alpha = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#FCAA67", alpha = 0.6) +  #the old pink: #edafb8
  geom_line(linewidth = 1.2, colour = "#FCAA67") +
  labs(
    title = expression("(a) "*italic("Daucus carota")*" "),
    x = "Hoverfly abundance",
    y = "Plant  reproductive success"
  ) +
  theme_classic() +
  theme(
    # Increase title and axis text size
    plot.title = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.ticks = element_line(color = "black")
  )
a

###########################
# (b) Radish: fruit set  #
###########################

#Show the effect for the social bees
newdat_mod.1_Radish <- data.frame(A_socialBees_Radish.dayly = seq(0, 31, by = 1))
newdat_mod.1_Radish$A_socialBees_Radish.dayly.z <- (newdat_mod.1_Radish$A_socialBees_Radish.dayly-mean(df7b$A_socialBees_Radish.dayly))/sd(df7b$A_socialBees_Radish.dayly)
newdat_mod.1_Radish$A_Apis_Radish.dayly.z <- 0 
newdat_mod.1_Radish$A_solitaryBees_Radish.dayly.z <- 0 
#newdat_mod.1_Radish$A_otherAculeata_Radish.dayly.z <- 0
newdat_mod.1_Radish$A_Syrphidae_Radish.dayly.z <- 0 
#newdat_mod.1_Radish$A_Coleoptera_Radish.dayly.t.z <-0

Xmat_mod.1_Radish <- model.matrix(~ A_Apis_Radish.dayly.z +  A_socialBees_Radish.dayly.z + A_solitaryBees_Radish.dayly.z +  A_Syrphidae_Radish.dayly.z , data=newdat_mod.1_Radish) #A_otherAculeata_Radish.dayly.z + A_Coleoptera_Radish.dayly.t.z
fitmat_mod.1_Radish <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Radish))
for(i in 1:nsim)fitmat_mod.1_Radish[,i] <- plogis(Xmat_mod.1_Radish %*% bsim.mod.1_Radish@fixef[i,])#plogis for binary models
newdat_mod.1_Radish$lwr <- apply(fitmat_mod.1_Radish, 1, quantile, prob=0.025)
newdat_mod.1_Radish$upr <- apply(fitmat_mod.1_Radish, 1, quantile, prob=0.975)
newdat_mod.1_Radish$fit <-plogis(Xmat_mod.1_Radish %*% fixef(mod.1_Radish))

df7b$Fruit_set <- df7b$n_flowers_with_fruits / (df7b$n_flowers_with_fruits + df7b$n_flowers_without_fruits)

b <- ggplot(newdat_mod.1_Radish, aes(x = A_socialBees_Radish.dayly, y = fit)) +
  geom_point(data = df7b, aes(x = A_socialBees_Radish.dayly, y = Fruit_set), 
                        color = "black", 
                        shape = 16, 
                        fill = "black", 
                        size = 4,
                        stroke = 0, alpha = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#D33665", alpha = 0.6) +
  geom_line(size = 1.2, colour = "#D33665") +
  labs(
    title = expression("(b) "*italic("Raphanus sativus")*" "),
    x = "Social bee abundance",
    y = "Plant  reproductive success"
  ) +
  theme_classic() +
  theme(
    # Increase title and axis text size
    plot.title = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.ticks = element_line(color = "black")
  )


b

############################
# (c) Sainfoin: fruit set  #
############################

#Show the effect for the all bees
newdat_mod.1b_Sainfoin <- data.frame(A_allBees_Sainfoin.dayly = seq(0, 3, by = 0.01))
newdat_mod.1b_Sainfoin$A_allBees_Sainfoin.dayly.z <- (newdat_mod.1b_Sainfoin$A_allBees_Sainfoin.dayly-mean(df7d$A_allBees_Sainfoin.dayly))/sd(df7d$A_allBees_Sainfoin.dayly)

Xmat_mod.1b_Sainfoin <- model.matrix(~ A_allBees_Sainfoin.dayly.z, data=newdat_mod.1b_Sainfoin)
fitmat_mod.1b_Sainfoin <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1b_Sainfoin))
for(i in 1:nsim)fitmat_mod.1b_Sainfoin[,i] <- plogis(Xmat_mod.1b_Sainfoin %*% bsim.mod.1b_Sainfoin@fixef[i,])#plogis for binary models
newdat_mod.1b_Sainfoin$lwr <- apply(fitmat_mod.1b_Sainfoin, 1, quantile, prob=0.025)
newdat_mod.1b_Sainfoin$upr <- apply(fitmat_mod.1b_Sainfoin, 1, quantile, prob=0.975)
newdat_mod.1b_Sainfoin$fit <-plogis(Xmat_mod.1b_Sainfoin %*% fixef(mod.1b_Sainfoin))


df7d$Fruit_set <- df7d$n_flowers_with_fruits / (df7d$n_flowers_with_fruits + df7d$n_flowers_without_fruits)

c <- ggplot(newdat_mod.1b_Sainfoin, aes(x = A_allBees_Sainfoin.dayly, y = fit)) +
  geom_point(data = df7d, aes(x = A_allBees_Sainfoin.dayly, y = Fruit_set), 
             color = "black", 
             shape = 16, 
             fill = "black", 
             size = 4,
             stroke = 0, alpha = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#E08DAC", alpha = 0.6) +
  geom_line(size = 1.2, colour = "#E08DAC") +
  labs(
    title = expression("(c) "*italic("Onobrychis viciifolia")*" "),
    x = "Total bee abundance",
    y = "Plant  reproductive success"
  ) +
  theme_classic() +
  theme(
    # Increase title and axis text size
    plot.title = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.ticks = element_line(color = "black")
  )
  
c

###########################
# (d) Comfrey: fruit set  #
###########################

#Show the effect for the all bees
newdat_mod.1_Comfrey <- data.frame(A_Bombus_Comfrey.dayly = seq(0, 7, by = 0.01))
newdat_mod.1_Comfrey$A_Bombus_Comfrey.dayly.z <- (newdat_mod.1_Comfrey$A_Bombus_Comfrey.dayly-mean(df7e$A_Bombus_Comfrey.dayly))/sd(df7e$A_Bombus_Comfrey.dayly)

Xmat_mod.1_Comfrey <- model.matrix(~ A_Bombus_Comfrey.dayly.z, data=newdat_mod.1_Comfrey)
fitmat_mod.1_Comfrey <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Comfrey))

for(i in 1:nsim)fitmat_mod.1_Comfrey[,i] <- plogis(Xmat_mod.1_Comfrey %*% bsim.mod.1_Comfrey@fixef[i,])#plogis for binary models

newdat_mod.1_Comfrey$lwr <- apply(fitmat_mod.1_Comfrey, 1, quantile, prob=0.025)
newdat_mod.1_Comfrey$upr <- apply(fitmat_mod.1_Comfrey, 1, quantile, prob=0.975)
newdat_mod.1_Comfrey$fit <-plogis(Xmat_mod.1_Comfrey %*% fixef(mod.1_Comfrey))

df7e$Fruit_set <- df7e$n_flowers_with_seeds / (df7e$n_flowers_with_seeds + df7e$n_flowers_without_seeds)

d <- ggplot(newdat_mod.1_Comfrey, aes(x = A_Bombus_Comfrey.dayly, y = fit)) +
  geom_point(data = df7e, aes(x = A_Bombus_Comfrey.dayly, y = Fruit_set), 
             color = "black", 
             shape = 16, 
             fill = "black", 
             size = 4,
             stroke = 0, alpha = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#8D80AD", alpha = 0.6) +
  geom_line(size = 1.2, colour = "#8D80AD") +
  labs(
    title = expression("(d) "*italic("Symphytum officinale")*" "),
    x = "Bumblebee abundance",
    y = "Plant  reproductive success"
  ) +
  theme_classic() +
  theme(
    # Increase title and axis text size
    plot.title = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.ticks = element_line(color = "black")
  )
d

# now the coefficient plots

# Load dataset
df_CoefPlot <- read.table("results/beta_coefficients_reproductive_success_all.txt", 
                          header = TRUE, sep = ";", dec = ".")

# Define effect order once
effect_order <- c("Honeybees", "social_Bees", "solitary_Bees", "all_Bees", "Wasps", "Hoverflies", "Beetles")
# Define order in sentence case for correct label factor
effect_order <- c("Honeybees", "social_Bees", "solitary_Bees", "all_Bees", "Wasps", "Hoverflies", "Beetles")

library(stringr)

effect_order_labels <- str_to_sentence(gsub("_", " ", effect_order))  # Now matches your label format

df_CoefPlot <- df_CoefPlot %>%
  mutate(
    fixeff = recode(fixeff,
                    "Bumblebees" = "social_Bees",
                    "solitary_Bees" = "solitary_Bees"),
    Sig = ifelse(lower > 0 | upper < 0, "Credible", "Uncertain"),
    label = factor(str_to_sentence(gsub("_", " ", fixeff)),
                   levels = rev(effect_order_labels)),
    shape = case_when(  # ← NEW
      Response == "Fruit_set" ~ 16,  # circle
      Response == "Seed_set" ~ 17    # triangle
    )
  )


df_carrot_forest <- df_CoefPlot %>%
  filter(Phytometer == "Carrot", Effect == "Abundance", Response == "Seed_set", fixeff != "Intercept") %>%
  mutate(
    fixeff = factor(fixeff, levels = rev(effect_order))  # only if needed for internal use
  )

forest_a <- ggplot(df_carrot_forest, aes(x = mean, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(aes(color = Sig), size = 4,   shape =   17) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = Sig), height = 0, linewidth = 0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  #coord_flip()+
  labs(
    x = "Mean parameter estimate\nand 95% credible interval",
    y = "",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black")
  )

forest_a


#radish
df_radish_combined <- df_CoefPlot %>%
  filter(Phytometer == "Radish",
         Effect == "Abundance",
         Response %in% c("Fruit_set", "Seed_set"),
         fixeff != "Intercept")

forest_b <- ggplot(df_radish_combined, aes(x = mean, y = label, group = Response)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  
  # ERROR BARS: must use group and position_dodge
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = Sig),
    height = 0,
    linewidth = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  
  # POINTS: also use group and matching dodge
  geom_point(
    aes(color = Sig, shape = Response),
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  
  scale_shape_manual(values = c("Fruit_set" = 16, "Seed_set" = 17)) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  
  labs(
    x = "Mean parameter estimate\nand 95% credible interval",
    y = "",
    shape = "Response",
    color = "Sig"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black")
  )

forest_b

#sain foin

df_sainfoin_combined <- df_CoefPlot %>%
  filter(
    Phytometer == "Sainfoin",
    Effect == "Abundance",
    Response %in% c("Fruit_set", "Seed_set"),
    fixeff != "Intercept"
  )

forest_c <- ggplot(df_sainfoin_combined, aes(x = mean, y = label, group = Response)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = Sig),
    height = 0,
    linewidth = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  
  geom_point(
    aes(color = Sig, shape = Response),
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  
  scale_shape_manual(values = c("Fruit_set" = 16, "Seed_set" = 17)) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  
  labs(
    x = "Mean parameter estimate\nand 95% credible interval",
    y = "",
    shape = "Response",
    color = "Sig"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black")
  )

forest_c

df_comfrey_combined <- df_CoefPlot %>%
  filter(
    Phytometer == "Comfrey",
    Effect == "Abundance",
    Response %in% c("Fruit_set", "Seed_set"),
    fixeff != "Intercept"
  )

forest_d <- ggplot(df_comfrey_combined, aes(x = mean, y = label, group = Response)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = Sig),
    height = 0,
    linewidth = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  
  geom_point(
    aes(color = Sig, shape = Response),
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  
  scale_shape_manual(values = c("Fruit_set" = 16, "Seed_set" = 17)) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  
  labs(
    x = "Mean parameter estimate\nand 95% credible interval",
    y = "",
    shape = "Response",
    color = "Sig"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.x = element_line(color = "black")
  )

forest_d

library(patchwork)

final_layout <- (
  ((a | forest_a) + plot_layout(widths = c(1.2, 0.7))) /
    ((b | forest_b) + plot_layout(widths = c(1.2, 0.7))) /
    ((c | forest_c) + plot_layout(widths = c(1.2, 0.7))) /
    ((d | forest_d) + plot_layout(widths = c(1.2, 0.7)))
)

final_layout #figure s5


# save as 20 x 10 inches portrait, rest of it cleaned up in inkscape

##Run until here##

############################################################################################################################
############################################################################################################################
############################################################################################################################

