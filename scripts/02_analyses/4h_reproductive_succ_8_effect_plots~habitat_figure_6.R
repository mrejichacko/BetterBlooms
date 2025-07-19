#####################################################################
###
### 4h. Pollinator experiment: 
### Modelling plant reproductive success: effects plots
### Code by: David Frey and Merin Reji Chacko
### Last edited: 18.07.2025
### 
#####################################################################
#####################################################################
#####################################################################
###
### Prep

###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
### load 3 a and 4a-g! it's not enough to just load 3a!!!!!!!!!!!!###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###

#########################
# (a) Carrot: seed set  #
#########################

#1. high plant species richness
newdat_mod.1_Carrot.1 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Carrot.1$PlantS.z <- max(df7a$PlantS.z) 
newdat_mod.1_Carrot.1$Urban_500.z <- (newdat_mod.1_Carrot.1$Urban_500-mean(df7a$Urban_500))/sd(df7a$Urban_500)
Xmat_mod.1_Carrot.1 <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Carrot.1)
newdat_mod.1_Carrot.1$fit <- exp(Xmat_mod.1_Carrot.1 %*% fixef(mod.1_Carrot))# exp = inverse link function for poisson models. 
fitmat_mod.1_Carrot.1 <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Carrot.1))
for(i in 1:nsim) fitmat_mod.1_Carrot.1[,i] <- exp(Xmat_mod.1_Carrot.1 %*% bsim.mod.1_Carrot@fixef[i,])
newdat_mod.1_Carrot.1$lwr <- apply(fitmat_mod.1_Carrot.1, 1, quantile, prob=0.025)
newdat_mod.1_Carrot.1$upr <- apply(fitmat_mod.1_Carrot.1, 1, quantile, prob=0.975)

#2. low plant species richness
newdat_mod.1_Carrot.2 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Carrot.2$PlantS.z <- min(df7a$PlantS.z) 
newdat_mod.1_Carrot.2$Urban_500.z <- (newdat_mod.1_Carrot.2$Urban_500-mean(df7a$Urban_500))/sd(df7a$Urban_500)

Xmat_mod.1_Carrot.2 <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Carrot.2)

newdat_mod.1_Carrot.2$fit <- exp(Xmat_mod.1_Carrot.2 %*% fixef(mod.1_Carrot))# exp = inverse link function for poisson models. 

fitmat_mod.1_Carrot.2 <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Carrot.2))
for(i in 1:nsim) fitmat_mod.1_Carrot.2[,i] <- exp(Xmat_mod.1_Carrot.2 %*% bsim.mod.1_Carrot@fixef[i,])
newdat_mod.1_Carrot.2$lwr <- apply(fitmat_mod.1_Carrot.2, 1, quantile, prob=0.025)
newdat_mod.1_Carrot.2$upr <- apply(fitmat_mod.1_Carrot.2, 1, quantile, prob=0.975)

library(ggplot2)

df7a$Seed_set <- df7a$n_seeds  # optional, just for clarity

a <- ggplot() +
  # Observed data points
  geom_point(data = df7a, aes(x = Urban_500, y = Seed_set), 
             shape = 17, color = "black", size = 4, alpha = 1) +
  
  # High richness ribbon + line
  geom_ribbon(data = newdat_mod.1_Carrot.1, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#1B7837", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Carrot.1, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#1B7837") +
  
  # Low richness ribbon + line
  geom_ribbon(data = newdat_mod.1_Carrot.2, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#A6611A", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Carrot.2, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#A6611A") +
  
  labs(
    title = expression("(a) "*italic("Daucus carota")*""),
    x = "Densification",
    y = "Pollination success"
  ) +
  scale_x_continuous(limits = c(0.19, 0.84)) +
  theme_classic(base_size=14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.line = element_line(linewidth = 0.8)
  )
a

##########################
# (b) Radish: fruit set  #
##########################

#1. high plant species richness
newdat_mod.1_Radish.1 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Radish.1$Urban_500.z <- (newdat_mod.1_Radish.1$Urban_500-mean(df7b$Urban_500))/sd(df7b$Urban_500)
newdat_mod.1_Radish.1$PlantS.z <- max(df7b$PlantS.z) 

Xmat_mod.1_Radish.1 <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Radish.1) #
fitmat_mod.1_Radish.1 <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Radish.1))
for(i in 1:nsim)fitmat_mod.1_Radish.1[,i] <- plogis(Xmat_mod.1_Radish.1 %*% bsim.mod.1_Radish@fixef[i,])#plogis for binary models
newdat_mod.1_Radish.1$lwr <- apply(fitmat_mod.1_Radish.1, 1, quantile, prob=0.025)
newdat_mod.1_Radish.1$upr <- apply(fitmat_mod.1_Radish.1, 1, quantile, prob=0.975)
newdat_mod.1_Radish.1$fit <-plogis(Xmat_mod.1_Radish.1 %*% fixef(mod.1_Radish))

#2. low plant species richness
newdat_mod.1_Radish.2 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Radish.2$Urban_500.z <- (newdat_mod.1_Radish.2$Urban_500-mean(df7b$Urban_500))/sd(df7b$Urban_500)
newdat_mod.1_Radish.2$PlantS.z <- min(df7b$PlantS.z) 

Xmat_mod.1_Radish.2 <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Radish.2) #
fitmat_mod.1_Radish.2 <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Radish.2))
for(i in 1:nsim)fitmat_mod.1_Radish.2[,i] <- plogis(Xmat_mod.1_Radish.2 %*% bsim.mod.1_Radish@fixef[i,])#plogis for binary models
newdat_mod.1_Radish.2$lwr <- apply(fitmat_mod.1_Radish.2, 1, quantile, prob=0.025)
newdat_mod.1_Radish.2$upr <- apply(fitmat_mod.1_Radish.2, 1, quantile, prob=0.975)
newdat_mod.1_Radish.2$fit <-plogis(Xmat_mod.1_Radish.2 %*% fixef(mod.1_Radish))

df7b$Fruit_set <- df7b$n_flowers_with_fruits / (df7b$n_flowers_with_fruits + df7b$n_flowers_without_fruits)

b <- ggplot() +
  geom_point(data = df7b, aes(x = Urban_500, y = Fruit_set), 
             shape = 16, color = "black", size = 4, alpha = 1) +
  
  geom_ribbon(data = newdat_mod.1_Radish.1, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#1B7837", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Radish.1, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#1B7837") +
  
  geom_ribbon(data = newdat_mod.1_Radish.2, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#A6611A", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Radish.2, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#A6611A") +
  
  labs(
    title = expression("(b) "*italic("Raphanus sativus")*""),
    x = "Densification",
    y = "Pollination success"
  ) +
  scale_x_continuous(limits = c(0.19, 0.84)) +
  theme_classic(base_size=14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8)
  )

b

###########################
# (d) Comfrey: fruit set  #
###########################

#1. high plant species richness
newdat_mod.1_Comfrey.1 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Comfrey.1$Urban_500.z <- (newdat_mod.1_Comfrey.1$Urban_500-mean(df7e$Urban_500))/sd(df7e$Urban_500)
newdat_mod.1_Comfrey.1$PlantS.z <- max(df7e$PlantS.z) 

Xmat_mod.1_Comfrey <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Comfrey.1)
fitmat_mod.1_Comfrey <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Comfrey.1))
for(i in 1:nsim)fitmat_mod.1_Comfrey[,i] <- plogis(Xmat_mod.1_Comfrey %*% bsim.mod.1_Comfrey@fixef[i,])#plogis for binary models
newdat_mod.1_Comfrey.1$lwr <- apply(fitmat_mod.1_Comfrey, 1, quantile, prob=0.025)
newdat_mod.1_Comfrey.1$upr <- apply(fitmat_mod.1_Comfrey, 1, quantile, prob=0.975)
newdat_mod.1_Comfrey.1$fit <-plogis(Xmat_mod.1_Comfrey %*% fixef(mod.1_Comfrey))

#2. high plant species richness
newdat_mod.1_Comfrey.2 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Comfrey.2$Urban_500.z <- (newdat_mod.1_Comfrey.2$Urban_500-mean(df7e$Urban_500))/sd(df7e$Urban_500)
newdat_mod.1_Comfrey.2$PlantS.z <- min(df7e$PlantS.z) 

Xmat_mod.1_Comfrey <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Comfrey.2)
fitmat_mod.1_Comfrey <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Comfrey.2))
for(i in 1:nsim)fitmat_mod.1_Comfrey[,i] <- plogis(Xmat_mod.1_Comfrey %*% bsim.mod.1_Comfrey@fixef[i,])#plogis for binary models
newdat_mod.1_Comfrey.2$lwr <- apply(fitmat_mod.1_Comfrey, 1, quantile, prob=0.025)
newdat_mod.1_Comfrey.2$upr <- apply(fitmat_mod.1_Comfrey, 1, quantile, prob=0.975)
newdat_mod.1_Comfrey.2$fit <-plogis(Xmat_mod.1_Comfrey %*% fixef(mod.1_Comfrey))

df7e$Fruit_set <- df7e$n_flowers_with_seeds / (df7e$n_flowers_with_seeds + df7e$n_flowers_without_seeds)

d <- ggplot() +
  geom_point(data = df7e, aes(x = Urban_500, y = Fruit_set), 
             shape = 16, color = "black", size = 4, alpha = 1) +
  
  geom_ribbon(data = newdat_mod.1_Comfrey.1, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#1B7837", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Comfrey.1, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#1B7837") +
  
  geom_ribbon(data = newdat_mod.1_Comfrey.2, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#A6611A", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Comfrey.2, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#A6611A") +
  
  labs(
    title = expression("(d) "*italic("Symphytum officinale")),
    x = "Densification",
    y = "Pollination success"
  ) +
  theme_classic(base_size=14) +
  scale_x_continuous(limits = c(0.19, 0.84)) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8)
  )

d


############################
# (c) Sainfoin: fruit set  #
############################

#1. high plant species richness
newdat_mod.1_Sainfoin.1 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Sainfoin.1$Urban_500.z <- (newdat_mod.1_Sainfoin.1$Urban_500-mean(df8$Urban_500))/sd(df8$Urban_500)
newdat_mod.1_Sainfoin.1$PlantS.z <- max(df8$PlantS.z) 

Xmat_mod.1_Sainfoin.1 <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Sainfoin.1)
fitmat_mod.1_Sainfoin.1 <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Sainfoin.1))
for(i in 1:nsim)fitmat_mod.1_Sainfoin.1[,i] <- plogis(Xmat_mod.1_Sainfoin.1 %*% bsim.mod.1_Sainfoin@fixef[i,])#plogis for binary models
newdat_mod.1_Sainfoin.1$lwr <- apply(fitmat_mod.1_Sainfoin.1, 1, quantile, prob=0.025)
newdat_mod.1_Sainfoin.1$upr <- apply(fitmat_mod.1_Sainfoin.1, 1, quantile, prob=0.975)
newdat_mod.1_Sainfoin.1$fit <-plogis(Xmat_mod.1_Sainfoin.1 %*% fixef(mod.1_Sainfoin))

#2. low plant species richness
newdat_mod.1_Sainfoin.2 <- data.frame(Urban_500 = seq(0, 1, by = 0.01))
newdat_mod.1_Sainfoin.2$Urban_500.z <- (newdat_mod.1_Sainfoin.2$Urban_500-mean(df8$Urban_500))/sd(df8$Urban_500)
newdat_mod.1_Sainfoin.2$PlantS.z <- min(df8$PlantS.z) 

Xmat_mod.1_Sainfoin.2 <- model.matrix(~ PlantS.z + Urban_500.z + PlantS.z:Urban_500.z, data=newdat_mod.1_Sainfoin.2)
fitmat_mod.1_Sainfoin.2 <-  matrix(ncol = nsim, nrow = nrow(newdat_mod.1_Sainfoin.2))
for(i in 1:nsim)fitmat_mod.1_Sainfoin.2[,i] <- plogis(Xmat_mod.1_Sainfoin.2 %*% bsim.mod.1_Sainfoin@fixef[i,])#plogis for binary models
newdat_mod.1_Sainfoin.2$lwr <- apply(fitmat_mod.1_Sainfoin.2, 1, quantile, prob=0.025)
newdat_mod.1_Sainfoin.2$upr <- apply(fitmat_mod.1_Sainfoin.2, 1, quantile, prob=0.975)
newdat_mod.1_Sainfoin.2$fit <-plogis(Xmat_mod.1_Sainfoin.2 %*% fixef(mod.1_Sainfoin))

df8$Fruit_set <- df8$n_flowers_with_fruits / (df8$n_flowers_with_fruits + df8$n_flowers_without_fruits)

c <- ggplot() +
  geom_point(data = df8, aes(x = Urban_500, y = Fruit_set), 
             shape = 16, color = "black", size = 4, alpha = 1) +
  
  geom_ribbon(data = newdat_mod.1_Sainfoin.1, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#1B7837", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Sainfoin.1, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#1B7837") +
  
  geom_ribbon(data = newdat_mod.1_Sainfoin.2, 
              aes(x = Urban_500, ymin = lwr, ymax = upr), 
              fill = "#A6611A", alpha = 0.3) +
  geom_line(data = newdat_mod.1_Sainfoin.2, 
            aes(x = Urban_500, y = fit), size = 1.2, color = "#A6611A") +
  
  labs(
    title = expression("(c) "*italic("Onobrychis viciifolia")),
    x = "Densification",
    y = "Pollination success"
  ) +
  theme_classic(base_size=14) +
  scale_x_continuous(limits = c(0.19, 0.84)) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8)
  )

c

library(grid)
library(patchwork)


a | b | c | d
#for now: save landscape 10.5 x 10

# adding the coefficient plots!

unique(coeff.mods.ReprSucc.GxL.table$Model)
coef_for_plot <- subset(coeff.mods.ReprSucc.GxL.table,
                        Model %in% c("mod.1_Carrot", "mod.1_Radish", "mod.1_Sainfoin", "mod.1_Comfrey"))
coef_for_plot

unique(coef_for_plot$fixeff)


# Optional: make fixeff an ordered factor so terms plot in the desired order
library(dplyr)
library(forcats)

coef_for_plot <- coef_for_plot %>%
  filter(fixeff != "Intercept") %>%
  mutate(
    fixeff = fct_recode(fixeff,
                        "Local flowering\nplant species richness" = "Garden",
                        "Landscape-scale\nhabitat loss" = "Urban_500.z",
                        "Interaction" = "PlantS.z:Urban_500.z"
    ),
    fixeff = factor(fixeff, levels = (c(
      "Local flowering\nplant species richness",
      "Landscape-scale\nhabitat loss",
      "Interaction"
    )))
  )


coef_for_plot

coef_for_plot <- coef_for_plot %>%
  mutate(
    Effect = ifelse(lower < 0 & upper > 0, "Uncertain", "Credible")
  )

a_coef <- ggplot(subset(coef_for_plot, Phytometer == "Carrot"), 
                 aes(x = mean, y = fixeff, color = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(shape = 17, size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  coord_flip() +
  labs(
    x = "Mean parameter estimate and 95% credible interval",
    y = "For gridding purpose - delete!",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(hjust = 0, color = "black", size = 16),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "italic", size = 14)
  )
a_coef

b_coef <- ggplot(subset(coef_for_plot, Phytometer == "Radish"), 
                 aes(x = mean, y = fixeff, color = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(shape = 16, size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  coord_flip() +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "italic", size = 14)
  )

b_coef

c_coef <- ggplot(subset(coef_for_plot, Phytometer == "Sainfoin"), 
                 aes(x = mean, y = fixeff, color = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(shape = 16, size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  coord_flip() +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "italic", size = 14)
  )
c_coef

d_coef <- ggplot(subset(coef_for_plot, Phytometer == "Comfrey"), 
                 aes(x = mean, y = fixeff, color = Effect)) +
  coord_flip()+
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E01A4F", linewidth = 0.8) +
  geom_point(shape = 16, size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.8) +
  scale_color_manual(values = c("Credible" = "black", "Uncertain" = "lightgrey")) +
  #coord_flip() +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(hjust = 0, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "italic", size = 14)
  )

d_coef

library(patchwork)

top_row <- a + a_coef + b + b_coef + 
  plot_layout(widths = c(2, 1, 2, 1))  # Main plots = 2, Coef plots = 1

bottom_row <- c + c_coef + d + d_coef + 
  plot_layout(widths = c(2, 1, 2, 1))

# Stack the rows
final_plot <- top_row / bottom_row
final_plot
