#Anova
#https://www.datanovia.com/en/lessons/ancova-in-r/
# ANCOVA makes several assumptions about the data, such as:
# 
# 1. Linearity between the covariate and the outcome variable at each level of the
# grouping variable.
# 
# This can be checked by creating a grouped scatter plot of the covariate and the outcome variable.
# 
# 2. Homogeneity of regression slopes. The slopes of the regression lines, formed
# by the covariate and the outcome variable, should be the same for each group.
# 
# This assumption evaluates that there is no interaction between the outcome and
# the covariate. The plotted regression lines by groups should be parallel.
# 
# 3. The outcome variable should be approximately normally distributed.
# This can be checked using the Shapiro-Wilk test of normality on the model residuals.
# 
# 4. Homoscedasticity or homogeneity of residuals variance for all groups.
# The residuals are assumed to have a constant variance (homoscedasticity)
# 5. No significant outliers in the groups


# one way ANOVA------

library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

# Data preparation----
# Load and prepare the data
data("anxiety", package = "datarium")
anxiety <- anxiety %>%
  select(id, group, t1, t3) %>%
  rename(pretest = t1, posttest = t3)
anxiety[14, "posttest"] <- 19
# Inspect the data by showing one random row by groups
set.seed(123)
anxiety %>% sample_n_by(group, size = 1)

# Check assumptions-----


# 1. Linearity assumption

ggscatter(
  anxiety, x = "pretest", y = "posttest",
  color = "group", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = group)
  )


# 2. Homogeneity of regression slopes

anxiety %>% anova_test(posttest ~ group*pretest)
#There was homogeneity of regression slopes as the interaction term was not
#statistically significant, F(2, 39) = 0.13, p = 0.88.


# 3. Normality of residuals

# Fit the model, the covariate goes first
model <- lm(posttest ~ pretest + group, data = anxiety)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid)


# if the Shapiro Wilk test was not significant (p > 0.05),  we can assume normality of residuals

# 4. Homogeneity of variances

model.metrics %>% levene_test(.resid ~ group)

# if the Levene’s test was not significant (p > 0.05),  we can assume
# homogeneity of the residual variances for all groups.

# 5. Outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# There were no outliers in the data, as assessed by no cases with standardized
# residuals greater than 3 in absolute value.


# Computation----

res.aov <- anxiety %>% anova_test(posttest ~ pretest + group)
get_anova_table(res.aov)

#After adjustment for pre-test anxiety score, there was a statistically
#significant difference in post-test anxiety score between the groups, F(2, 41)
#= 218.63, p < 0.0001.


# Post-hoc test----
# Pairwise comparisons can be performed to identify which groups are different.
# The Bonferroni multiple testing correction is applied. This can be easily done
# using the function emmeans_test() [rstatix package], a wrapper around the
# emmeans package, which needs to be installed. Emmeans stands for estimated
# marginal means (aka least square means or adjusted means).

# Pairwise comparisons
library(emmeans)
pwc <- anxiety %>% 
  emmeans_test(
    posttest ~ group, covariate = pretest,
    p.adjust.method = "bonferroni"
  )
pwc


# Display the adjusted means of each group
# Also called as the estimated marginal means (emmeans)
get_emmeans(pwc)




#Report----

# Visualization: line plots with p-values
pwc <- pwc %>% add_xy_position(x = "group", fun = "mean_se")
ggline(get_emmeans(pwc), x = "group", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(pwc, hide.ns = TRUE, tip.length = FALSE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )