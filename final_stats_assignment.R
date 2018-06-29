# final_stats_assignment.R
# Final statistical analyses for metabolic rate of Manduca sexta data
# for basic stats assignment due 29 June 2018
# Taryn Joy Jacobs, Jessica Kipling, Carlin Landsberg
# 18 June 2018

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readr)
library(ggpubr)
library(e1071)
library(corrplot)
library(RColorBrewer)
library(grid)

# Load data ---------------------------------------------------------------

ms_rate <- read_delim("ms_rate.csv", ";", 
                      escape_double = FALSE, trim_ws = TRUE)

# imported textfile from desktop
View(ms_rate)

# Having a look at the data -----------------------------------------------

# summarise
summary(ms_rate)

# Density plot ------------------------------------------------------------

# Visualise the denisty of the data
ggplot(data = ms_rate, aes(x = LogMrate, fill = as.factor(Instar))) + 
  geom_density(alpha = 0.4)
  # data seems like its not skewed for each instar 

# Skewness & Kurtosis -----------------------------------------------------

# skewness, kurtosis...
ms_rate %>% 
  group_by(Instar) %>% 
  summarise(mean_wt = mean(BodySize),
            median_wt = median(BodySize),
            skew_wt = skewness(BodySize),
            kurt_wt = kurtosis(BodySize)) 
  # 1 and 4 kurtosis

# Boxplot visualisation ---------------------------------------------------

# Boxplot
ggplot(data = ms_rate, aes(x = as.factor(Instar), y = Mrate)) +
  geom_boxplot(aes(fill = as.factor(Instar)), notch = TRUE)
  # Mostly looks significant (notches not overlapping, clearly for 3 and 4. 1 and 2 a bit unclear)

# Data analyses (first checking assumptions) -----------------------------------------------------------------

# check assumptions instar and rate

# normality
ms_rate %>% 
  group_by(Instar) %>% 
  summarise(r_norm_dist = as.numeric(shapiro.test(LogMrate)[2])) 
  # All data normal, except for instsar 5. Removed instar 5 from raw data?

# homoscedasticity 
ms_rate %>% 
  group_by(Instar) %>% 
  summarise(r_norm_dist = as.numeric(shapiro.test(LogMrate)[2]),
            r_norm_var = var(LogMrate)) 
 
ms_rate %>% 
  summarise(body_norm_dist = as.numeric(shapiro.test(BodySize)[2]),
            body_norm_var = var(BodySize))
  # body size not normally distributed

ms_rate %>% 
  summarise(body_norm_dist = as.numeric(shapiro.test(LogBodySize)[2]),
            body_norm_var = var(LogBodySize))
  # log body size also not normally distributed

# Note, we cannot perform a t-test becasue more than two factors (4 instars)
# so we can try an ANOVA ... single-factor anova

# ANOVA Metabolic rate ----------------------------------------------------

# H0: Metabolic rate is not higher for larger instars
# H1: Metabolic rate is higher for larger instars

# single factor ANOVA for instar and metabolic rate
mr.aov <- aov(LogMrate ~ as.factor(Instar), data = ms_rate)
summary(mr.aov)
  # ----- RESULTS -----
#                   Df Sum Sq Mean Sq F value Pr(>F)    
# as.factor(Instar)   3  76.40  25.467   177.5 <2e-16 ***
#  Residuals         252  36.16   0.143                   
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1   
  # ----- ----- -----
    # Pr < 0.05 so we reject the null hypothesis, therefore
      # Metabolic rate is higher for larger instars

# Visualise using a boxplot
mr_box <- ggplot(data = ms_rate, aes(x = as.factor(Instar), y = LogMrate, fill = as.factor(Instar))) +
  geom_boxplot(notch = TRUE) +
  annotate("text", x = 1, y = 0.1, label = "a") +
  annotate("text", x = 2, y = 0.85, label = "b") +
  annotate("text", x = 3, y = 1.5, label = "c") +
  annotate("text", x = 4, y = 1.83, label = "d") +
  labs(x = "Instar", y = "Metabolic rate", fill = "Instar",
       caption = "Fig. 1. Notched boxplots showing the relationship between instar stage and metabolic rate of the tobacco hornworm caterpillar, Manduca sexta.") +
  theme(plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0)) +
  theme_classic()
mr_box

   # to be used in results section
ggsave(plot = mr_box, filename = "Metabolic boxplot final2.png")
# notches do not overlap so must be significantly different?
# this is just a visualisation, so lets check it with a statistical analysis

# Post-hoc test
TukeyHSD(mr.aov)
plot(TukeyHSD(mr.aov))
  # none cross the zero mark... 
  # p-adj < 0.05 (p-adj = 0 for all interactions) for all so does this mean they're all significant?
  # ie all instars have significantly different metabolic rates?

# Correlation -------------------------------------------------------------

# subset dataset because we are looking at only instar, body size and metabolic rate
met_sub <- ms_rate %>% 
  select(-X1, -Computer, -BodySize, -CO2ppm, -Mrate)

# H0: There is no relationship between body size and metabolic rate
# H1: There is a relationship between body size and metabolic rate

# visualise the data (scatterplot)
plot(ms_rate$LogBodySize, ms_rate$LogMrate)

# normality for log body size and log metabolic rate
met_norm <- met_sub %>% 
  gather(key = "variable") %>% 
  group_by(variable) %>% 
  summarise(variable_norm = as.numeric(shapiro.test(value)[2]))
met_norm
  # but data is not narmal so we use Kendall

# Kendall rank correlation
corr <- cor(met_sub, method = "kendall")
corr

# visualise results form these correlations
    # to be used in results section
plot_panel <- (corrplot(corr, type = "upper",
                        tl.col = "black", addCoef.col = "grey70"))
  # overall high correlatation coefficients for all variables, therefore a strong, positive relationship exists

# correlation between body size and metabolic rate (Kendall)
cor.test(ms_rate$LogBodySize, ms_rate$LogMrate, method = "kendall")
# p < 0.05, so we reject the null hypothesis, therefore
# there is a relationship between body size and metabolic rate
# tau = 0.8158754, so about 82%% of this relationship is explained

# visualise plot
tau_print <- paste0("tau = ", 
                  round(cor(x = met_sub$LogBodySize, met_sub$LogMrate, method = "kendall"),2))

plot_corr <- ggplot(data = ms_rate, aes(x = LogBodySize, y = LogMrate)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = -2.2, y = 1.7, label = tau_print) +
  labs(x = "Body mass (g)", y = "Metabolic rate", caption = "Fig. 2. Correlation showing the relationship between body mass (g) and metabolic rate of Manduca sexta. 
") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5)) +
  theme_classic()
plot_corr

    # to be used in results section
ggsave(plot = plot_corr, filename = "correlation plot final2.png")

# Regression --------------------------------------------------------------

# H0: Metabolic rate is not influenced by body size
# H1: Metabolic rate is influenced by body size 

met_lm <- lm(LogBodySize ~ LogMrate, data = ms_rate)
summary(met_lm)
  # pr < 0.05 so we reject the null hypothesis, therefore
  # Metabolic rate is influenced by body size
  # adjusted r-sq value is 0.9235 showing that approx 92% of variance is explained

slope <- round(met_lm$coef[2], 3)
# p.val <- round(coefficients(summary(met_lm))[2, 4], 3) # approx 0 so...
p.val <- 0.001
r2 <- round(summary(met_lm)$r.squared, 3)
# ggplot(data = trans_met_co2, aes(x = Instar, y = LogMrate)) +
#  geom_point() +
#  stat_smooth(method = "lm")

plot_lm <- ggplot(data = ms_rate, aes(x = LogBodySize, y = LogMrate)) +
  geom_point() +
  annotate("text", x = -2.5, y = 1.9, label = paste0("slope == ", slope, "~(min/min)"), parse = TRUE, hjust = 0) +
  annotate("text", x = -2.5, y = 1.7, label = paste0("italic(p) < ", p.val), parse = TRUE, hjust = 0) +
  annotate("text", x = -2.5, y = 1.5, label = paste0("italic(r)^2 == ", r2), parse = TRUE, hjust = 0) +
  stat_smooth(method = "lm") +
  labs(x = "Body mass (g)", y = "Metabolic rate", caption = "Fig. 3. Linear regression showing the influence of body mass (g) on metabolic rate of Manduca sexta.
") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5)) +
  theme_classic()
plot_lm

    # to be used in results section
ggsave(plot = plot_lm, filename = "regression plot final.png")

# regression for QPAS... see papers..
ms_rate_scale <- ms_rate %>% 
  mutate(scale = (Mrate)^(-1/4))
ms_rate_scale

scale_lm <- lm(LogBodySize ~ scale, data = ms_rate_scale)
summary(scale_lm)
  # pr < 0.05 so we reject the null hypothesis, therefore
  # QPAS (Metabolic rate^ -1/4) is influenced by body size
  # adjusted r-sq value is 0.8597 showing that approx 86% of variance is explained

slope1 <- round(scale_lm$coef[2], 3)
# p.val <- round(coefficients(summary(met_lm))[2, 4], 3) # approx 0 so...
p.val <- 0.001
r2.1 <- round(summary(scale_lm)$r.squared, 3)

scale_plot_lm <- ggplot(data = ms_rate_scale, aes(x = LogBodySize, y = scale)) +
  geom_point() +
  annotate("text", x = -0.7, y = 2.2, label = paste0("slope == ", slope1, "~(min/min)"), parse = TRUE, hjust = 0) +
  annotate("text", x = -0.7, y = 2.0, label = paste0("italic(p) < ", p.val), parse = TRUE, hjust = 0) +
  annotate("text", x = -0.7, y = 1.8, label = paste0("italic(r)^2 == ", r2.1), parse = TRUE, hjust = 0) +
  stat_smooth(method = "lm") +
  labs(x = "Body mass (g)", y = "Metabolic rate (scale)", caption = "Fig. 5. Linear regression of quarter-power allometric scaling (QPAS) relationship of body mass and metabolic rate
") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5)) +
  theme_classic()
scale_plot_lm

ggsave(plot = scale_plot_lm, filename = "scale plot regression final.png")
