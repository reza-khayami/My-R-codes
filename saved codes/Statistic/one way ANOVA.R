#Anova
#https://www.datanovia.com/en/lessons/anova-in-r/
# The ANOVA test makes the following assumptions about the data:
#
# 1. Independence of the observations.
# Each subject should belong to only one group. There is no relationship between
# the observations in each group. Having repeated measures for the same participants is not allowed.

# 2. No significant outliers in any cell of the design 

# 3. Normality. the data for each design cell should be approximately normally distributed. 

# 4. Homogeneity of variances. The variance of the outcome variable should be
# equal in every cell of the design. Before computing ANOVA test, you need to
# perform some preliminary tests to check if the assumptions are met.

# Note that, if the above assumptions are not met there are a non-parametric alternative
# (Kruskal-Wallis test) to the one-way ANOVA.

# one way ANOVA------

library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

# Data preparation----
# Load and prepare the data
data("PlantGrowth")
set.seed(1234)
PlantGrowth %>% sample_n_by(group, size = 1)
levels(PlantGrowth$group)
PlantGrowth <- PlantGrowth %>%
  reorder_levels(group, order = c("ctrl", "trt1", "trt2"))

# The one-way ANOVA can be used to determine whether the means plant growths are
# significantly different between the three conditions.


# Summary statistics-----

PlantGrowth %>%
  group_by(group) %>%
  get_summary_stats(weight, type = "mean_sd")

# Visualization
ggboxplot(PlantGrowth, x = "group", y = "weight")


# Check assumptions----

#1.Outliers
PlantGrowth %>% 
  group_by(group) %>%
  identify_outliers(weight)

# Note that, in the situation where you have extreme outliers, this can be due
# to: 1) data entry errors, measurement errors or unusual values.
#
# Yo can include the outlier in the analysis anyway if you do not believe the
# result will be substantially affected. This can be evaluated by comparing the
# result of the ANOVA test with and without the outlier.
#
# It’s also possible to keep the outliers in the data and perform robust ANOVA
# test using the WRS2 package.


# 2. Normality assumption



# The normality assumption can be checked by using one of the following two
# approaches:
#
# 1. Analyzing the ANOVA model residuals to check the normality for all groups
# together. This approach is easier and it’s very handy when you have many
# groups or if there are few data points per group.

# 2. Check normality for each group separately. This approach might be used when
# you have only a few groups and many data points per group. In this section,
# we’ll show you how to proceed for both option 1 and 2.



# 2.1 Check normality assumption by analyzing the model residuals.

# QQ plot and Shapiro-Wilk test of normality are used. QQ plot draws the
# correlation between a given data and the normal distribution.

# Build the linear model
model  <- lm(weight ~ group, data = PlantGrowth)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


# In the QQ plot, as all the points fall approximately along the reference line,
# we can assume normality. This conclusion is supported by the Shapiro-Wilk
# test. If the p-value is not significant (here p = 0.13),  we can assume normality.

# 2.2 Check normality assumption by groups.
# Computing Shapiro-Wilk test for each group level. If the data is normally
# distributed, the p-value should be greater than 0.05.


PlantGrowth %>%
  group_by(group) %>%
  shapiro_test(weight)


# Note that, if your sample size is greater than 50, the normal QQ plot is
# preferred because at larger sample sizes the Shapiro-Wilk test becomes very
# sensitive even to a minor deviation from normality.
#
# QQ plot draws the correlation between a given data and the normal
# distribution. Create QQ plots for each group level:
#
ggqqplot(PlantGrowth, "weight", facet.by = "group")


# 3. Homogneity of variance assumption
plot(model, 1)

# In the plot above, there is no evident relationships between residuals and
# fitted values (the mean of each groups), which is good. So, we can assume the
# homogeneity of variances.


# It’s also possible to use the Levene’s test to check the homogeneity of variances:
  PlantGrowth %>% levene_test(weight ~ group)

  
  # From the output above, we can see that the p-value is > 0.05, which is not
  # significant. This means that, there is not significant difference between
  # variances across groups. Therefore, we can assume the homogeneity of
  # variances in the different treatment groups.
  #
  # In a situation where the homogeneity of variance assumption is not met, you
  # can compute the Welch one-way ANOVA test using the function
  # welch_anova_test()[rstatix package]. This test does not require the
  # assumption of equal variances.
  
  # Computation ----
  res.aov <- PlantGrowth %>% anova_test(weight ~ group)
  res.aov
  
  # From the above ANOVA table, it can be seen that there are significant
  # differences between groups (p = 0.016), which are highlighted with “*“, F(2,
  # 27) = 4.85, p = 0.16, eta2[g] = 0.26.
  
  
  # Post-hoc tests---- A significant one-way ANOVA is generally followed up by
  # Tukey post-hoc tests to perform multiple pairwise comparisons between
  # groups. Key R function: tukey_hsd() [rstatix].
  
  # Pairwise comparisons
  pwc <- PlantGrowth %>% tukey_hsd(weight ~ group)
  pwc
  
  
  # The output contains the following columns:
  #   
  # estimate: estimate of the difference between means of the two groups
  # conf.low, conf.high: the lower and the upper end point of the confidence interval at 95% (default)
  # p.adj: p-value after adjustment for the multiple comparisons.
  
  
  # Report ----
  
  # Visualization: box plots with p-values
  pwc <- pwc %>% add_xy_position(x = "group")
  ggboxplot(PlantGrowth, x = "group", y = "weight") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  
  
  
  # Relaxing the homogeneity of variance assumption ----
  
  # The classical one-way ANOVA test requires an assumption of equal variances
  # for all groups. In our example, the homogeneity of variance assumption
  # turned out to be fine: the Levene test is not significant.
  #
  # How do we save our ANOVA test, in a situation where the homogeneity of
  # variance assumption is violated?
  #
  # The Welch one-way test is an alternative to the standard one-way ANOVA in
  # the situation where the homogeneity of variance can’t be assumed (i.e.,
  # Levene test is significant). In this case, the Games-Howell post hoc test or
  # pairwise t-tests (with no assumption of equal variances) can be used to
  # compare all possible combinations of group differences.
  
  
  # Welch One way ANOVA test
  res.aov2 <- PlantGrowth %>% welch_anova_test(weight ~ group)
  
  # Pairwise comparisons (Games-Howell)
  pwc2 <- PlantGrowth %>% games_howell_test(weight ~ group)
  
  
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "group", step.increase = 1)
  ggboxplot(PlantGrowth, x = "group", y = "weight") +
    stat_pvalue_manual(pwc2, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(res.aov2, detailed = TRUE),
      caption = get_pwc_label(pwc2)
    )
  
  # You can also perform pairwise comparisons using pairwise t-test with no
  # assumption of equal variances:
    
    pwc3 <- PlantGrowth %>% 
    pairwise_t_test(
      weight ~ group, pool.sd = FALSE,
      p.adjust.method = "bonferroni"
    )
  pwc3
  