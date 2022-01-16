#ggplot
#if outcome is not 0,1
data2 <- data %>% 
  mutate(prob = ifelse(Status == "Cancer", 1, 0))

mtcars %>% 
  ggplot(aes(hp, vs)) +
  geom_point(aes(color = factor(vs))) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("Logistic regression model fit") +
  xlab("BiomarkerS") 


#base r
fit = glm(vs ~ hp, data=mtcars, family=binomial)
newdat <- data.frame(hp=seq(min(mtcars$hp), max(mtcars$hp),len=100))
newdat$vs = predict(fit, newdata=newdat, type="response")
plot(vs~hp, data=mtcars, col="red4")
lines(vs ~ hp, newdat, col="green4", lwd=2)