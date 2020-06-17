## SYSTEM ##
## Packages used in this script ##############
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "splines", "brms", "readxl", "mgcv", 
               "merTools", "gridExtra", "grid", "ggplot2", "microbenchmark",
               "car", "lme4", "lmtest", "foreach", "doParallel","gvlma")
set.seed(112)
############# Loading data ###################
setwd("C:/Users/Eier/Downloads")
patients_crp_org <- read_excel('crp_pasient.xlsx')
controls_crp_org <- read_excel('crp_control.xlsx')
##############################################
### Parallel cores used ######################
library(doParallel)
clusterobj<-makeCluster(4)
registerDoParallel(clusterobj)
##############################################

## Produces sample space of newdata needed in predictions ##
get.newdata <- function(method.A, method.B)
{
  seq(from = min(method.B), to = max(method.B), by = (max(method.B) - min(method.A))*0.01)
}

## REAL DATA - USED IN CHAPTER 2.1 ##
#a######## Data: Architect vs Advia ###########################
patients.arc.adv <- dplyr::select(patients_crp_org, sample, replicat, "Architect", "Advia") %>%
  na.omit %>% rename(A = "Architect", B = "Advia") %>%
  mutate_at(vars(sample), as.factor) %>%
  mutate(mm = 0.5*(A+B), ld = log(A) - log(B))

controls.arc.adv <- dplyr::select(controls_crp_org, sample, replicat, "Architect", "Advia") %>% 
  na.omit %>% rename(A = "Architect", B = "Advia") %>%
  mutate_at(vars(sample), as.factor) %>%
  mutate(mm = 0.5*(A+B), ld = log(A) - log(B))

## MOR ##
patients.arc.adv.MOR <- patients.arc.adv %>% group_by(sample) %>% summarise(mA = mean(A), mB = mean(B))
controls.arc.adv.MOR <- controls.arc.adv %>% group_by(sample) %>% summarise(mA = mean(A), mB = mean(B))

#b######## Data: Dimension vs Cobas #####################################
patients.dim.cob <- dplyr::select(patients_crp_org, sample, replicat, "Dimension", "Cobas") %>%
  na.omit %>% rename(A = "Dimension", B = "Cobas") %>%
  mutate_at(vars(sample), as.factor) %>%
  mutate(mm = 0.5*(A+B), ld = log(A) - log(B))
  
controls.dim.cob <- dplyr::select(controls_crp_org, sample, replicat, "Dimension", "Cobas") %>% 
  na.omit %>% rename(A = "Dimension", B = "Cobas") %>%
  mutate_at(vars(sample), as.factor) %>%
  mutate(mm = 0.5*(A+B), ld = log(A) - log(B))

## MOR ##
patients.dim.cob.MOR <- patients.dim.cob %>% group_by(sample) %>% summarise(mA = mean(A), mB = mean(B))
controls.dim.cob.MOR <- controls.dim.cob %>% group_by(sample) %>% summarise(mA = mean(A), mB = mean(B))


#c######## Data: Advia vs. Dimension ############################################
patients.adv.dim <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Advia", "Dimension") %>% # Select relevant components of main df
  na.omit %>% rename(A = "Advia", B = "Dimension") %>% # rename methods as A and B for convenience
  mutate_at(vars(sample), as.factor) %>% # From numeric to factor so that we can use that later
  mutate(mm = (A + B)/2, ld = log(A) - log(B)) # Creating two new coloumns using BA-transformation

controls.adv.dim <- dplyr::select(controls_crp_org, sample, replicat, "Advia", "Dimension") %>% 
  na.omit %>% rename(A = "Advia", B = "Dimension") %>%
  mutate_at(vars(sample), as.factor) %>%
  mutate(mm = (A + B)/2, ld = log(A) - log(B))

## MOR ##
patients.adv.dim.MOR <- patients.adv.dim %>% group_by(sample) %>% summarise(mA = mean(A), mB = mean(B))
controls.adv.dim.MOR <- controls.adv.dim %>% group_by(sample) %>% summarise(mA = mean(A), mB = mean(B))
## the MOR data frames are the means of the replicates within samples ##
 
#2##### Methods of evaluation 1 : Ordinary least squares regression (OLSR) #########
## Linear model objects (all replicates) ##
olsr.arc.adv <- lm(data = patients.arc.adv, formula = A ~ B)
olsr.dim.cob <- lm(data = patients.dim.cob, formula = A ~ B)
olsr.adv.dim <- lm(data = patients.adv.dim, formula = A ~ B)

## Linear model objects (mean of replicates) ##
olsr.arc.adv.MOR <- lm(data = patients.arc.adv.MOR, formula = mA ~ mB)
olsr.dim.cob.MOR <- lm(data = patients.dim.cob.MOR, formula = mA ~ mB)
olsr.adv.dim.MOR <- lm(data = patients.adv.dim.MOR, formula = mA ~ mB)
se.coef(olsr.arc.adv.MOR)
se.coef(olsr.arc.adv)
## New data ##
newdata.aa <- get.newdata(patients.arc.adv$A,patients.arc.adv$B)
newdata.dc <- get.newdata(patients.dim.cob$A,patients.dim.cob$B)
newdata.ad <- get.newdata(patients.adv.dim$A,patients.adv.dim$B)

## Prediction intervals - All replicates ##
pred.arc.adv <- data.frame(new = newdata.aa, predict(object = olsr.arc.adv, newdata = list(B = newdata.aa), interval = "prediction", level = 0.99))
pred.dim.cob <- data.frame(new = newdata.dc, predict(object = olsr.dim.cob, newdata = list(B = newdata.dc), interval = "prediction", level = 0.99))
pred.adv.dim <- data.frame(new = newdata.ad, predict(object = olsr.adv.dim, newdata = list(B = newdata.ad), interval = "prediction", level = 0.99))

## Prediction intervals - Mean of replicates ##
pred.arc.adv.MOR <- data.frame(new = newdata.aa, predict(object = olsr.arc.adv.MOR, newdata = list(mB = newdata.aa), interval = "prediction", level = 0.99))
pred.dim.cob.MOR <- data.frame(new = newdata.dc, predict(object = olsr.dim.cob.MOR, newdata = list(mB = newdata.dc), interval = "prediction", level = 0.99))
pred.adv.dim.MOR <- data.frame(new = newdata.ad, predict(object = olsr.adv.dim.MOR, newdata = list(mB = newdata.ad), interval = "prediction", level = 0.99))

#3##### Plots with 99% prediction bands ################################

plot2a <- ggplot() +
  geom_ribbon(data = pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  geom_line(data = pred.arc.adv, aes(x = new, y = fit), size = 2, color = "violet", alpha = 0.5) +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color  = "blue",size = 2) +
  geom_point(data = controls.arc.adv, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  labs(title = "Architect vs. Advia", subtitle = "With 99% prediction bands") +
  ylab("Architect") + xlab("Advia")

# Dimension vs Cobas
plot2b <- ggplot() +
  geom_ribbon(data = pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  geom_line(data = pred.dim.cob, aes(x = new, y = fit), size = 2, color = "violet", alpha = 0.5) +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), size = 2, color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  labs(title = "Dimension vs. Cobas", subtitle = "With 99% prediction bands") +
  ylab("Dimension") + xlab("Cobas")

# Advia vs dimension
plot2c <- ggplot() +
  geom_ribbon(data = pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", size = 1, alpha = 0.3, color = "black") +
  geom_line(data = pred.adv.dim, aes(x = new, y = fit), color = "violet", size = 2, alpha = 0.5) +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  labs(title = "Advia vs. Dimension", subtitle = "With 99% prediction bands") +
  ylab("Advia") +
  xlab("Dimension")

plot(plot2a)
plot(plot2b)
plot(plot2c)

grid.arrange(plot2a,plot2b,plot2c)


plot2d <- ggplot() +
  geom_ribbon(data = pred.arc.adv.MOR, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  geom_line(data = pred.arc.adv.MOR, aes(x = new, y = fit), color = "violet", alpha = 0.5, size = 2) +
  geom_point(data = patients.arc.adv.MOR, aes(x = mB, y = mA), color  = "blue", size = 2) +
  geom_point(data = controls.arc.adv.MOR, aes(x = mB, y = mA, shape = sample), size = 3, color = "red") +
  labs(title = "OLSR with MOR", subtitle = "With 99% prediciton bands") +
  ylab("Architect") + xlab("Advia")

# Dimension vs Cobas
plot2e <- ggplot() +
  geom_ribbon(data = pred.dim.cob.MOR, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  #geom_ribbon(data = pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "orange", alpha = 0.2, color = "black") +
  geom_line(data = pred.dim.cob.MOR, aes(x = new, y = fit), color = "violet", size = 2, alpha = 0.5) +
  geom_point(data = patients.dim.cob.MOR, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.dim.cob.MOR, aes(x = mB, y = mA, shape = sample), size = 3, color = "red") +
  labs(title = "OLSR with MOR", subtitle = "With 99% prediciton bands") +
  ylab("Dimension") + xlab("Cobas")

# Advia vs dimension
plot2f <- ggplot() +
  geom_ribbon(data = pred.adv.dim.MOR, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  #geom_ribbon(data = pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "orange", alpha = 0.2, color = "black") +
  geom_line(data = pred.adv.dim.MOR, aes(x = new, y = fit), color = "violet", size = 2, alpha = 0.5) +
  geom_point(data = patients.adv.dim.MOR, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.adv.dim.MOR, aes(x = mB, y = mA, shape = sample), size = 3, color = "red") +
  labs(title = "OLSR with MOR", subtitle = "Green ribbon is 99% prediciton bands") +
  ylab("Advia") + xlab("Dimension")

plot(plot2d)
plot(plot2e)
plot(plot2f)

grid.arrange(plot2d,plot2e,plot2f)

#4##### Residual plots for plot OLSR and formal tests ##########################
plot2a.res <-  ggplot(data = NULL, mapping = aes(x = fitted(olsr.arc.adv), y = residuals(olsr.arc.adv))) +
  geom_point() + geom_hline(yintercept = 0) +
  ylab("Residuals") + xlab("Fitted") + 
  labs(title = "Architect vs. Advia", subtitle = "Residuals vs. fitted plot")
plot2b.res <-  ggplot(data = NULL, mapping = aes(x = fitted(olsr.dim.cob), y = residuals(olsr.dim.cob))) +
  geom_point() + geom_hline(yintercept = 0) +
  ylab("Residuals") + xlab("Fitted") + 
  labs(title = "Dimension vs. Cobas", subtitle = "Residuals vs. fitted plot")
plot2c.res <-  ggplot(data = NULL, mapping = aes(x = fitted(olsr.adv.dim), y = residuals(olsr.adv.dim))) +
  geom_point() + geom_hline(yintercept = 0) +
  ylab("Residuals") + xlab("Fitted") + 
  labs(title = "Advia vs. Dimension", subtitle = "Residuals vs. fitted plot")
plot2d.res <-  ggplot(data = NULL, mapping = aes(x = fitted(mean.olsr.arc.adv), y = residuals(mean.olsr.arc.adv))) +
  geom_point() + geom_hline(yintercept = 0) +
  ylab("Residuals") + xlab("Fitted") + 
  labs(title = "Architect vs. Advia", subtitle = "Residuals vs. fitted plot (MOR)")
plot2e.res <-  ggplot(data = NULL, mapping = aes(x = fitted(mean.olsr.dim.cob), y = residuals(mean.olsr.dim.cob))) +
  geom_point() + geom_hline(yintercept = 0) +
  ylab("Residuals") + xlab("Fitted") + 
  labs(title = "Dimension vs. Cobas", subtitle = "Residuals vs. fitted plot (MOR)")
plot2f.res <-  ggplot(data = NULL, mapping = aes(x = fitted(mean.olsr.adv.dim), y = residuals(mean.olsr.adv.dim))) +
  geom_point() + geom_hline(yintercept = 0) +
  ylab("Residuals") + xlab("Fitted") + 
  labs(title = "Advia vs. Dimension", subtitle = "Residuals vs. fitted plot (MOR)")
res.plots <- grid.arrange(plot2a.res, plot2b.res, plot2c.res, plot2d.res, plot2e.res, plot2f.res, nrow = 2, ncol = 3)

# Normality
shapiro.test(residuals(olsr.arc.adv))
shapiro.test(residuals(olsr.dim.cob))
shapiro.test(residuals(olsr.adv.dim))
shapiro.test(residuals(mean.olsr.arc.adv))
shapiro.test(residuals(mean.olsr.dim.cob))
shapiro.test(residuals(mean.olsr.adv.dim))

# Homoscedasticy
bptest(olsr.arc.adv)
bptest(olsr.dim.cob)
bptest(olsr.adv.dim)
bptest(mean.olsr.arc.adv)
bptest(mean.olsr.dim.cob)
bptest(mean.olsr.adv.dim)

# Autocorrelation
dwtest(olsr.arc.adv)
dwtest(olsr.dim.cob)
dwtest(olsr.adv.dim)
dwtest(mean.olsr.arc.adv)
dwtest(mean.olsr.dim.cob)
dwtest(mean.olsr.adv.dim)

#5##### Log-log plots (linearity?)##########################################################

# Architect vs Advia
plot3a <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "Log-log plot") + ylab("log(Architect)") + xlab("log(Advia)")

# Dimension vs Cobas
plot3b <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "Log-log plot") + ylab("ln(Dimension)") + xlab("ln(Cobas)")

# Advia vs Dimension
plot3c <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "Log-log plot") + ylab("ln(Advia)") + xlab("ln(Dimension)")

plot(plot3a); plot(plot3b); plot(plot3c)

#6##### Log-log plots with olsr and prediction bands ############################

log.olsr.arc.adv <- lm(data = patients.arc.adv, formula = log(A) ~ log(B))
log.olsr.dim.cob <- lm(data = patients.dim.cob, formula = log(A) ~ log(B))
log.olsr.adv.dim <- lm(data = patients.adv.dim, formula = log(A) ~ log(B))

# Prediction intervals
log.pred.arc.adv <- data.frame(new = log(newdata.aa), predict(object = log.olsr.arc.adv, newdata = list(B = newdata.aa), interval = "prediction", level = 0.99)) 
log.pred.dim.cob <- data.frame(new = log(newdata.dc), predict(object = log.olsr.dim.cob, newdata = list(B = newdata.dc), interval = "prediction", level = 0.99))
log.pred.adv.dim <- data.frame(new = log(newdata.ad), predict(object = log.olsr.adv.dim, newdata = list(B = newdata.ad), interval = "prediction", level = 0.99))

# Architect vs. Advia
plot4a <- ggplot() +
  geom_ribbon(data = log.pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  geom_line(data = log.pred.arc.adv, aes(x = new, y = fit), color = "violet", size = 2, alpha = 0.5) +
  geom_point(data = patients.arc.adv, aes(x = log(B), y = log(A)), color  = "blue", size = 2) +
  geom_point(data = controls.arc.adv, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 3) +
  labs(title = "ln(Architect) vs. ln(Advia)", subtitle = "With 99% prediction bands") +
  ylab("ln(Architect)") + xlab("ln(Advia)")

# Dimension vs Cobas
plot4b <- ggplot() +
  geom_ribbon(data = log.pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  geom_line(data = log.pred.dim.cob, aes(x = new, y = fit), color = "violet", alpha = 0.5, size = 2) +
  geom_point(data = patients.dim.cob, aes(x = log(B), y = log(A)), color  = "blue", size = 2) +
  geom_point(data = controls.dim.cob, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 3) +
  labs(title = "ln(Dimension) vs. ln(Cobas)", subtitle = "With 99% prediction bands") +
  ylab("ln(Dimension)") + xlab("ln(Cobas)")

# Advia vs dimension
plot4c <- ggplot() +
  geom_ribbon(data = log.pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 1) +
  geom_line(data = log.pred.adv.dim, aes(x = new, y = fit), color = "violet", size = 2, alpha = 0.5) +
  geom_point(data = patients.adv.dim, aes(x = log(B), y = log(A)), color  = "blue", size = 2) +
  geom_point(data = controls.adv.dim, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 3) +
  labs(title = "ln(Advia) vs. ln(Dimension)", subtitle = "With 99% prediction bands") +
  ylab("ln(Advia)") + xlab("ln(Dimension)")


# Plots
plot(plot4a)
plot(plot4b)
plot(plot4c)


grid.arrange(plot4a,plot2a, nrow = 1)
grid.arrange(plot4b,plot2b,plot4c,plot2c, nrow = 2)

# Comparison between ordinary plots and log-log plots
grid.arrange(plot2a,plot4a, nrow = 1)
grid.arrange(plot2b,plot4b, nrow = 1)
grid.arrange(plot2c,plot4c, nrow = 1)
#7##### Residual plots and formal tests ########################

plot5a <- ggplot(mapping = aes(x = fitted.values(olsr.arc.adv), y = residuals(olsr.arc.adv))) +
  geom_point(color = "black", shape = 17, size = 2) + geom_hline(yintercept = 0, size = 1) +
  labs(title = "Architect vs. Advia", subtitle = "Ordinary residual plot") +
  ylab("Residuals") + xlab("Fitted values")

plot5b <- ggplot(mapping = aes(x = fitted.values(olsr.dim.cob), y = residuals(olsr.dim.cob)), color = "black", shape = 18) +
  geom_point(color = "black", shape = 17, size = 2) +geom_hline(yintercept = 0, size = 1) +
  labs(title = "Dimension vs. Cobas", subtitle = "Ordinary residual plot") +
  ylab("Residuals") + xlab("Fitted values")

plot5c <- ggplot(mapping = aes(x = fitted.values(olsr.adv.dim), y = residuals(olsr.adv.dim)), color = "black", shape = 18) +
  geom_point(color = "black", shape = 17, size = 2) + geom_hline(yintercept = 0, size = 1) +
  labs(title = "Advia vs. Dimension", subtitle = "Ordinary residual plot") +
  ylab("Residuals") + xlab("Fitted values")

plot5d <- ggplot(mapping = aes(x = fitted.values(log.olsr.arc.adv), y = residuals(log.olsr.arc.adv)), color = "black") +
  geom_point() + geom_hline(yintercept = 0, size = 1) +
  labs(title = "Architect vs. Advia", subtitle = "log-log residual plot") +
  ylab("Residuals") + xlab("Fitted values")

plot5e <- ggplot(mapping = aes(x = fitted.values(log.olsr.dim.cob), y = residuals(log.olsr.dim.cob)), color = "black") +
  geom_point() + geom_hline(yintercept = 0, size = 1) +
  labs(title = "Dimension vs. Cobas", subtitle = "log-log residual plot") +
  ylab("Residuals") + xlab("Fitted values")

plot5f <- ggplot(mapping = aes(x = fitted.values(log.olsr.adv.dim), y = residuals(log.olsr.adv.dim)), color = "black") +
  geom_point() + geom_hline(yintercept = 0, size = 1) +
  labs(title = "Advia vs. Dimension", subtitle = "log-log residual plot") +
  ylab("Residuals") + xlab("Fitted values")

plot(plot5a); plot(plot5b); plot(plot5c)
plot(plot5d); plot(plot5e); plot(plot5f)
grid.arrange(plot5a, plot5d, nrow = 1)
grid.arrange(plot5b, plot5e, nrow = 1)
grid.arrange(plot5c, plot5f, nrow = 1)

##### Formal tests
# Normality
s1<-shapiro.test(residuals(log.olsr.arc.adv))
s2<-shapiro.test(residuals(log.olsr.dim.cob))
s3<-shapiro.test(residuals(log.olsr.adv.dim))

# Homoscedasticy
b1<-bptest(log.olsr.arc.adv)
b2<-bptest(log.olsr.dim.cob)
b3<-bptest(log.olsr.adv.dim)

# Autocorrelation
d1<-dwtest(log.olsr.arc.adv)
d2<-dwtest(log.olsr.dim.cob)
d3<-dwtest(log.olsr.adv.dim)


#7##### Methods of evaluation 2 :  Bland-Altman and polynomial regression (Needs work)

## Fixed effects models ##
ba.arc.adv <- lm(data = patients.arc.adv, formula = ld ~ poly(mm,4))
ba.dim.cob <- lm(data = patients.dim.cob, formula = ld ~ poly(mm,4))
ba.adv.dim <- lm(data = patients.adv.dim, formula = ld ~ poly(mm,4))

## Newdata ##
ba.new.aa <- get.newdata(patients.arc.adv$mm,patients.arc.adv$mm)
ba.new.dc <- get.newdata(patients.dim.cob$mm,patients.dim.cob$mm)
ba.new.ad <- get.newdata(patients.adv.dim$mm,patients.adv.dim$mm)

## Random effects models ##
rba.arc.adv <- lmer(data = patients.arc.adv, formula = ld ~ poly(mm,4) + (1|replicat))
rba.dim.cob <- lmer(data = patients.dim.cob, formula = ld ~ poly(mm,4) + (1|replicat))
rba.adv.dim <- lmer(data = patients.adv.dim, formula = ld ~ poly(mm,4) + (1|replicat))

## PREDICTIONS ##

## AA ##
rpred.ba.arc.adv <- data.frame(new = ba.new.aa, predictInterval(merMod = rba.arc.adv, newdata = data.frame(mm = ba.new.aa, replicat = 2), level = 0.99, n.sims = 9999))
pred.ba.arc.adv <- data.frame(new = ba.new.aa, predict(object = ba.arc.adv, newdata = list(mm = ba.new.aa), level = 0.99, interval = "prediction"))
## DC ##
rpred.ba.dim.cob <- data.frame(new = ba.new.dc, predictInterval(merMod = rba.dim.cob, newdata = data.frame(mm = ba.new.dc, replicat = 2), level = 0.99, n.sims = 9999))
pred.ba.dim.cob <- data.frame(new = ba.new.dc, predict(object = ba.dim.cob, newdata = list(mm = ba.new.dc), level = 0.99, interval = "prediction"))
## AD ##
rpred.ba.adv.dim <- data.frame(new = ba.new.ad, predictInterval(merMod = rba.adv.dim, newdata = data.frame(mm = ba.new.ad, replicat = 2), level = 0.99, n.sims = 9999, .parallel = T))
pred.ba.adv.dim <- data.frame(new = ba.new.ad, predict(object = ba.adv.dim, newdata = list(mm = ba.new.ad), level = 0.99, interval = "prediction"))

## PLOTS ##

## AA ##
plot.ba.aa <- ggplot() + geom_ribbon(data = rpred.ba.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.2, size = 1, color = "black") +
  geom_ribbon(data = pred.ba.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", size = 0.2) +
  geom_line(data = rpred.ba.arc.adv, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_line(data = rpred.ba.arc.adv, aes(x = new, y = fit), color = "gray", size = 1, linetype = 2) +
  geom_point(data = patients.arc.adv, aes(x = mm, y = ld), color = "blue", size = 2) +
  geom_point(data = controls.arc.adv, aes(x = mm, y = ld, shape = sample), color = "red", size = 3) +
  xlab("Mean of measurement methods") + ylab("Logarithmic difference between methods") +
  labs(title = "Bland-Altman assessment - AD - AR", subtitle = "Mixed effects model included")

## DC ##
plot.ba.dc <- ggplot() + geom_ribbon(data = rpred.ba.dim.cob, aes(x = new, ymin = lwr, ymax = upr), color = "black", fill = "green", size = 1, alpha = 0.2) +
  geom_ribbon(data = pred.ba.dim.cob, aes(x = new, ymin = lwr, ymax = upr), color = "black", fill = "green", size = 0.2, alpha = 0.3) +
  geom_line(data = rpred.ba.dim.cob, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_line(data = pred.ba.dim.cob, aes(x = new, y = fit), color = "gray", size = 1, linetype = 2) +
  geom_point(data = patients.dim.cob, aes(x = mm, y = ld), color = "blue", size = 2) +
  geom_point(data = controls.dim.cob, aes(x = mm, y = ld, shape = sample), color = "red", size = 3) +
  xlab("Mean of measurement methods") + ylab("Logarithmic difference between methods") +
  labs(title = "Bland-Altman assessment - DC - AR", subtitle = "Mixed effects model included")

## AD ##
plot.ba.ad <- ggplot() + geom_ribbon(data = rpred.ba.adv.dim, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "green", size = 1, alpha = 0.2) +
  geom_ribbon(data = pred.ba.adv.dim, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "green", size = 0.2, alpha = 0.3) +
  geom_line(data = rpred.ba.adv.dim, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_line(data = pred.ba.adv.dim, aes(x = new, y = fit), color = "gray", size = 1, linetype = 2) +
  geom_point(data = patients.adv.dim, aes(x = mm, y = ld), color = "blue", size = 2) +
  geom_point(data = controls.adv.dim, aes(x = mm, y = ld, shape = sample), color = "red", size = 3) +
  xlab("Mean of measurement methods") + ylab("Logarithmic difference between methods") +
  labs(title = "Bland-Altman assessment - AD - AR", subtitle = "Mixed effects model included")

plot(plot.ba.ad)

grid.arrange(plot.ba.aa,plot.ba.dc,plot.ba.ad)

#8##### Deming regression function ############################################
estimate.lambda <- function(method.A, method.B, replicates = 3)
{
  r <- replicates; n <- length(method.A) / r
  replicat <- sort(rep(1:r, n)); sample <- rep(1:n, r)
  df <- data.table::data.table(sample = sample, replicat = replicat, A = method.A, B = method.B) %>%
    group_by(sample) %>%
    mutate(mA = mean(A), mB = mean(B)) %>%
    mutate(YmM = (A - mA)^2, XmM = (B - mB)^2)
  sigma.ee <- (sum(df$YmM)) / (n*r - n)
  sigma.uu <- (sum(df$XmM)) / (n*r - n)
  return(lambda = sigma.ee / sigma.uu)
}

deming_fit <- function(method.A, method.B, R)
{
  A<-method.A;B<-method.B;lambda<-estimate.lambda(A,B,R);n<-length(A)
  print(A);print(B)
  S_AA <- crossprod(A-mean(A))[,]
  S_BB <- crossprod(B-mean(B))[,]
  S_BA <- crossprod(B-mean(B),A-mean(A))[,]

  b1 <- (S_AA-lambda * S_BB + sqrt((S_AA-lambda*S_BB)^2 + 4*lambda*S_BA^2))/(2*S_BA)
  b0 <- mean(A) - b1 * mean(B)
  
  sigma.uu <- ((S_AA + lambda*S_BB) - sqrt((S_AA-lambda*S_BB)^2 + 4*lambda*S_BA^2))/(2*lambda)/(n-1)
  s.vv <- crossprod(A-mean(A)-b1*(B-mean(B)))[,]/(n-2)
  var.b1 <- (S_BB*S_AA-S_BA^2)/(n*(S_BA/b1)^2)
  var.b0 <- s.vv/n + (mean(B))^2*var.b1
  #print(sigma.uu);print(s.vv)
  #print(var.b1);print(var.b0)
  return(list(lambda = lambda, coefficients = data.frame(b0=b0,b1=b1),
              sigma = sqrt(sigma.uu*(n-1)/(n-2)),
              V = rbind( c(var.b0, -(mean(B))*var.b1), c(-(mean(B))*var.b1, var.b1)),
              residuals = A - (b0+b1*B), fitted = b0+b1*B))
}

## Tests ##########################################
deming_fit(patients.adv.dim$A,patients.adv.dim$B,3)
deming_fit(patients.dim.cob$A,patients.dim.cob$B,3)
deming_fit(patients.arc.adv$A,patients.arc.adv$B,3)
###################################################

## PREDICTION BANDS - Parametric procedure ##

## FUNCTION - PREDICTION INTERVAL FOR ONE VALUE ##
deming.predict <-  function(method.A, method.B, newdata, level=0.99){
  n <- length(method.A)
  fit <- deming_fit(method.A,method.B,3)
  lambda <- fit$lambda
  sigma <- fit$sigma
  sigma.uu <- sigma^2*(n-2)/(n-1)
  V <- fit$V
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  ynew <- a+b*newdata
  xnew <- as.matrix(c(1,newdata))
  sigma.ee <- lambda*sigma.uu 
  t <- qt(1-(1-level)/2, n-2)
  ##  predict observed ynew given theoretical xnew
  sd.ynew.Fuller <- sqrt(sigma.ee + t(xnew)%*%V%*%xnew + (b^2+V[2,2])*sigma.uu)
  Lynew.Fuller <- ynew - t*sd.ynew.Fuller 
  Uynew.Fuller <- ynew + t*sd.ynew.Fuller 
  out <- data.table::data.table(new = newdata, fitted = ynew, lwr = Lynew.Fuller, upr = Uynew.Fuller) %>%
    rename(fit = "fitted.b0", lwr = "lwr.b0", upr = "upr.b0")
  out
}

## TEST OF FUNCTION - IS THE PREDICTED INTERVAL PLAUSABLE? ##
deming.predict(patients.adv.dim$A,patients.adv.dim$B,6)
deming.predict(patients.adv.dim$A,patients.adv.dim$B,9.33)
deming.predict(patients.arc.adv$A,patients.arc.adv$B,-33)
deming.predict(patients.dim.cob$A,patients.dim.cob$B,0)
#############################################################

## PREDICTION INTERVALS FOR MANY VALUES ##
deming.predictInterval <- function(method.A, method.B, newdata)
{
  foreach(i=1:length(newdata), .combine = rbind, .export = c("estimate.lambda","deming.predict","deming_fit"), .packages = "dplyr") %dopar% deming.predict(method.A = method.A, method.B = method.B, newdata = newdata[i])
}

## TEST OF FUNCTION - DO WE GET AN DATA-FRAME CONSISTING OF ALL INTERVALS? ##
## IS ALL PREDICTION INTERVALS PLAUSABLE? ##
## CAN WE USE BOTH INTEGERS AND FLOATS ? ##
deming.predictInterval(patients.adv.dim$A,patients.adv.dim$B,1:100)
deming.predictInterval(patients.dim.cob$A,patients.dim.cob$B, newdata = runif(101,5,100))


#9##### Methods of evaluation 3 : Deming regression prediction function #################################

##### Standard error with jackknife technique #####
get.se <- function(method.A,method.B,replicates)
{
  s.oo <- t(foreach(i=1:length(method.A), .combine = cbind, .export = c("deming.lm","estimate.lambda"), .packages = "dplyr") %dopar% t(deming.lm(method.A[-i],method.B[-i], replicates)$coefficients))
  s.org <- deming.lm(method.A,method.B,replicates)$coefficients; N<-length(method.A)
  s.org <- data.table::data.table(b0 = rep(s.org$b0,N),b1=rep(s.org$b1,N))  
  s.pse <- N * s.org - (N-1)*s.oo
  jfv <- sapply(s.pse,var)
  jfse <- sqrt(jfv/N)
  return(jfse)
}

## Deming regression assessment plots ##
dr.arc.adv <- deming.predictInterval(patients.arc.adv$A,patients.arc.adv$B,get.newdata(patients.arc.adv$A,patients.arc.adv$B))
dr.dim.cob <- deming.predictInterval(patients.dim.cob$A,patients.dim.cob$B,get.newdata(patients.dim.cob$A,patients.dim.cob$B))
dr.adv.dim <- deming.predictInterval(patients.adv.dim$A,patients.adv.dim$B,get.newdata(patients.adv.dim$A,patients.adv.dim$B))

plot.dr.aa <- ggplot() +
  geom_ribbon(data = dr.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", size = 1, alpha = 0.3, color = "black") +
  geom_line(data = dr.arc.adv, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_point(data = patients.arc.adv.MOR, aes(x = mB, y = mA), color = "blue", size = 2) +
  geom_point(data = controls.arc.adv.MOR, aes(x = mB, y = mA, shape = sample), color = "red", size = 3) +
  xlab("Advia") + ylab("Architect")+
  labs(title = "DR - AA - MOR", subtitle = "The green region is the 99% prediction bands")

plot.dr.dc <- ggplot() +
  geom_ribbon(data = dr.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", size = 1, alpha = 0.3, color = "black") +
  geom_line(data = dr.dim.cob, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_point(data = patients.dim.cob.MOR, aes(x = mB, y = mA), color = "blue", size = 2) +
  geom_point(data = controls.dim.cob.MOR, aes(x = mB, y = mA, shape = sample), color = "red", size = 3) +
  xlab("Cobas") + ylab("Dimension") +
  labs(title = "DR - DC - MOR", subtitle = "The green region is the 99% prediction bands")

plot.dr.ad <- ggplot() +
  geom_ribbon(data = dr.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", size = 1, alpha = 0.3, color = "black") +
  geom_line(data = dr.adv.dim, aes(x = new, y =  fit), color = "violet", size = 2) +
  geom_point(data = patients.adv.dim.MOR, aes(x = mB, y = mA), color = "blue", size = 2) +
  geom_point(data = controls.adv.dim.MOR, aes(x = mB, y = mA, shape = sample), color = "red", size = 3) +
  xlab("Dimension") + ylab("Advia") +
  labs(title = "DR- AD - MOR", subtitle = "The green region is the 99% prediction bands")

grid.arrange(plot.dr.aa,plot.dr.dc,plot.dr.ad)

dr.lm.aa <- deming_fit(patients.arc.adv$A,patients.arc.adv$B,3)
dr.lm.dc <- deming_fit(patients.dim.cob$A,patients.dim.cob$B,3)
dr.lm.ad <- deming_fit(patients.adv.dim$A,patients.adv.dim$B,3)

coef(dr.lm.ad)

## Residual plots - DR ##
dr.res.aa <- ggplot() + geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted(dr.lm.aa), y = residuals(dr.lm.aa)), size = 2) +
  ylab("Residuals") + xlab("Fitted values") +
  labs(title = "Residual plot - DR - AR", subtitle = "Architect vs. Advia")
dr.res.dc <- ggplot() + geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted(dr.lm.dc), y = residuals(dr.lm.dc)), size = 2) +
  ylab("Residuals") + xlab("Fitted values") +
  labs(title = "Residual plot - DR - AR", subtitle = "Dimension vs. Cobas")
dr.res.ad <- ggplot() + geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted(dr.lm.ad), y = residuals(dr.lm.ad)), size = 2) +
  ylab("Residuals") + xlab("Fitted values") +
  labs(title = "Residual plot - DR - AR", subtitle = "Advia vs. Dimension")

grid.arrange(dr.res.aa,dr.res.dc,dr.res.ad,ncol=2,nrow=2)













###################################
## Must be made manually I think
## Shapiro for DR ##

## BP for DR  ##

## DW for DR ##

#################################


##### PI with bootstrap technique #####
get.leverage <- function(method.B = NULL)
{
  design.matrix <- as.matrix(data.table::data.table(ones = rep(1, length(method.B)), B = method.B))  
  hat.matrix <- design.matrix%*%solve(t(design.matrix)%*%design.matrix)%*%t(design.matrix)
  leverage <- diag(hat.matrix)
  return(leverage)
}

get.sample <- function(method.A, method.B, replicates = 3, new)
{
  fit <- deming.lm(method.A = method.A, method.B = method.B, replicates = replicates) 
  y.p <- fit$coefficients$b0 + fit$coefficients$b1 * new
  residuals <- fit$residuals
  adjusted.residuals <- (residuals)/(sqrt(1 - get.leverage(method.B)))
  s <- adjusted.residuals - mean(adjusted.residuals)
  return(list(s = s, y.p = y.p, fit = fit))
}

bootstrap.resample <- function(s, fit, replicates = 3, new)
{
  s.bs <- sample(s,length(s),replace = T)
  y.bs.fit <- unlist(fitted(fit))+s.bs; r <- replicates
  x <- fit$model.frame$predictor
  fit.new <- deming.lm(y.bs.fit, x, replicates = r)
  residuals.new <- fit.new$residuals
  adjusted.residuals.new <- (residuals.new)/(sqrt(1 - get.leverage(x)))
  s.bs <- adjusted.residuals.new - mean(adjusted.residuals.new) 
  error.fit <- fit$coefficients$b0 - fit.new$coefficients$b0
  error.fit <- error.fit + (fit$coefficients$b1 - fit.new$coefficients$b1)*new
  return((unname(error.fit + sample(s.bs, size=1))))
}

bootstrap.pred <- function(method.A, method.B, resamples = 50, new)
{
  sample <- get.sample(method.A, method.B, replicates = 3, new = new)
  B <- resamples
  draws <- replicate(n = B, expr = bootstrap.resample(s = sample$s, fit = sample$fit, new = new))
  return(draws)
}

bootstrap.predictInterval <- function(method.A, method.B, resamples, level, upr, lwr)
{
  s <- get.sample(method.A,method.B,3,min(method.A))
  fit <- s$fit; s <- s$s; alpha<- 1-level;N<-resamples
  A<-method.A;B<-method.B; new<-lwr:upr
  df<- data.table::data.table(foreach(i=lwr:upr, .packages = c("dplyr", "tidyverse"), .export = c("bootstrap.pred","get.sample", "deming.lm","estimate.lambda","get.leverage", "bootstrap.resample"), .combine=cbind) %dopar% bootstrap.pred(A,B,N,i))
  df<- data.table::data.table(t(sapply(X=df,FUN=quantile,probs=c(alpha/2, 1 - alpha/2)))) %>%
    mutate(pred = fit$coefficients$b0 + fit$coefficients$b1 * new)  
  return(data.table::data.table(new = new, pred = df[,3], lwr=df[,1] + df[,3],upr=df[,2] + df[,3]))
}


## Newdata ##
get.newdata(patients.arc.adv$A,patients.arc.adv$B)
get.newdata(patients.dim.cob$A,patients.dim.cob$B)
get.newdata(patients.adv.dim$A,patients.adv.dim$B)

## Bootstrap tech ##
pidr.dim.cob <- bootstrap.predictInterval(patients.dim.cob$A,patients.dim.cob$B, 10^4, 0.99, 87.3, 1.5)
pidr.arc.adv <- bootstrap.predictInterval(patients.arc.adv$A,patients.arc.adv$B, 10^4, 0.99, 89.2, 1.4)
pidr.adv.dim <- bootstrap.predictInterval(patients.adv.dim$A,patients.adv.dim$B, 10^4, 0.99, 89.3, 1.5)

## Parametric tech ##
.pidr.dim.cob <- deming.predictInterval(patients.dim.cob$A,patients.dim.cob$B,newdata = seq(from = 1.5, to = 87.3, by = 0.5))
.pidr.arc.adv <- deming.predictInterval(patients.arc.adv$A,patients.arc.adv$B,newdata = seq(from = 1.4, to = 89.2, by = 0.5))
.pidr.adv.dim <- deming.predictInterval(patients.adv.dim$A,patients.adv.dim$B,newdata = seq(from = 1.5, to = 89.3, by = 0.5))


#10#### Plots of deming with 99% prediction interval ######################
plot6a <- ggplot() +
  geom_ribbon(data = .pidr.dim.cob, aes(x = new, ymin = lwr.V1, ymax = upr.V1), size = 1, fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = .pidr.dim.cob, aes(x=new,y=fitted), size = 2, color = "violet", alpha = 0.5) +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color = "blue", size = 2) +
  geom_point(data = controls.dim.cob, aes(x = B, y = A, shape=sample), size=3, color = "red") +
  xlab("Cobas") +
  ylab("Dimension") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

plot6b <- ggplot() +
  geom_ribbon(data = .pidr.arc.adv, aes(x = new, ymin = lwr.V1, ymax = upr.V1), size=1,fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = .pidr.arc.adv, aes(x = new, y = fitted), size = 2, color = "violet", alpha = 0.5) +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color = "blue", size = 2) +
  geom_point(data = controls.arc.adv, aes(x = B, y = A, shape=sample), size = 3, color = "red") +
  xlab("Advia") + ylab("Architect") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

plot6c <- ggplot() +
  geom_ribbon(data = .pidr.adv.dim, aes(x = new, ymin = lwr.V1, ymax = upr.V1), size=1, fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = .pidr.adv.dim, aes(x = new, y = fitted), size = 2, color = "violet", alpha = 0.5) +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color = "blue", size = 2) +
  geom_point(data = controls.adv.dim, aes(x = B, y = A,shape=sample), size = 3, color = "red") +
  xlab("Dimension") +
  ylab("Advia") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

# Plots - Deming regression + prediction bands

plot(plot6a)
plot(plot6b)
plot(plot6c)

grid.arrange(plot6a,plot6b)

###### Regression splines ######

# Models
spline.arc.adv.lm <- lm(data = patients.arc.adv, formula = A ~ ns(B, knots = c(30, 60)))
spline.dim.cob.lm <- lm(data = patients.dim.cob, formula = A ~ ns(B, knots = c(30, 60)))
spline.adv.dim.lm <- lm(data = patients.adv.dim, formula = A ~ ns(B, knots = c(30, 60)))
# predictions
spline.arc.adv.pred <- data.table::data.table(new = 3:95, predict(spline.arc.adv.lm,newdata = list(B=3:95),level = 0.99,interval = "prediction"))
spline.dim.cob.pred <- data.table::data.table(new = 3:95, predict(spline.dim.cob.lm,newdata = list(B=3:95),level = 0.99,interval = "prediction"))
spline.adv.dim.pred <- data.table::data.table(new = 3:95, predict(spline.adv.dim.lm,newdata = list(B=3:95),level = 0.99,interval = "prediction"))

# Plots
spline.plot.arc.adv <- ggplot() +
  geom_ribbon(data = spline.arc.adv.pred, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", size = 1, alpha = 0.3) +
  geom_line(data = spline.arc.adv.pred, aes(x = new, y = fit), color = "gray", alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(30, 60), size = 1, linetype = "dashed") +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color = "blue") +
  geom_point(data = controls.arc.adv, aes(x = B, y = A, shape = sample), color = "red") +
  xlab("Advia") + ylab("Architect") +
  labs(title = "Regression splines assessment evaluation",
       subtitle = "Architect vs. Advia (dashed lines are knots)")
spline.plot.dim.cob <- ggplot() +
  geom_ribbon(data = spline.dim.cob.pred, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", size = 1, alpha = 0.3) +
  geom_line(data = spline.dim.cob.pred, aes(x = new, y = fit), color = "gray", alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(30, 60), size = 1, linetype = "dashed") +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color = "blue") +
  geom_point(data = controls.dim.cob, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  xlab("Cobas") + ylab("Dimension") +
  labs(title = "Regression splines assessment evaluation",
       subtitle = "Dimension vs. Cobas (dashed lines are knots)")
spline.plot.adv.dim <- ggplot() +
  geom_ribbon(data = spline.adv.dim.pred, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", size = 1, alpha = 0.3) +
  geom_line(data = spline.adv.dim.pred, aes(x = new, y = fit), color = "gray", alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(30, 60), size = 1, linetype = "dashed") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color = "blue") +
  geom_point(data = controls.adv.dim, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  xlab("Dimension") + ylab("Advia") +
  labs(title = "Regression splines assessment evaluation",
       subtitle = "Advia vs. Dimension (dashed lines are knots)")

plot(spline.plot.arc.adv)
plot(spline.plot.dim.cob)
plot(spline.plot.adv.dim)

grid.arrange(spline.plot.dim.cob,spline.plot.adv.dim,
             nrow = 2, ncol = 1)





### Residual plots ###
res.plot.RS.arc.adv <- ggplot() +
  geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted.values(spline.arc.adv.lm), y = residuals(spline.arc.adv.lm))) +
  xlab("Fitted values") + ylab("Residuals") +
  labs(title = "Residual plot - RS", subtitle = "Architect vs. Advia")
res.plot.RS.dim.cob <- ggplot() +
  geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted.values(spline.dim.cob.lm), y = residuals(spline.dim.cob.lm))) +
  xlab("Fitted values") + ylab("Residuals") +
  labs(title = "Residual plot - RS", subtitle = "Dimension vs. Cobas")
res.plot.RS.adv.dim <- ggplot() +
  geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted.values(spline.adv.dim.lm), y = residuals(spline.adv.dim.lm))) +
  xlab("Fitted values") + ylab("Residuals") +
  labs(title = "Residual plot - RS", subtitle = "Advia vs. Dimension")

plot(res.plot.RS.arc.adv)
plot(res.plot.RS.dim.cob)
plot(res.plot.RS.adv.dim)

grid.arrange(res.plot.RS.arc.adv,res.plot.RS.dim.cob,res.plot.RS.adv.dim, nrow=1)

# Formal tests ##############################
shapiro.test(spline.arc.adv.lm$residuals)
shapiro.test(spline.dim.cob.lm$residuals)
shapiro.test(spline.adv.dim.lm$residuals)
#
bptest(spline.arc.adv.lm)
bptest(spline.dim.cob.lm)
bptest(spline.adv.dim.lm)
#
dwtest(spline.arc.adv.lm)
dwtest(spline.dim.cob.lm)
dwtest(spline.adv.dim.lm)





#############################################


#11#### Simulation ############################################

### We are simulating data ###
# can be used both for Patient Samples and Controls, and we can change param between them to test models
sim.data<- function(pairs, replicates, a=0, b=1.05, c=1, CVX, CVY, lower.limit, upper.limit)
{
  r <- replicates; n <- pairs
  sample <- (rep(1:n, each = r)) # Samples
  replicat <- rep(1:r, each = 1, times = n) # Replicates
  y.true <- rep(runif(n, lower.limit, upper.limit), each = r) # The sample space
  tmp <- data.table::data.table(sample, replicat, y.true) %>% 
    mutate(x.true = a * y.true ^ 2 + b * y.true + c) %>% # Some linear or second degree polynomial relationship may be added
    rowwise() %>% # We are using rows as aritmethics
    mutate(A = y.true * (1 + rnorm(1, 0, CVY))) %>%
    mutate(B = x.true * (1 + rnorm(1, 0, CVX))) %>%
    mutate(ld = log(A) - log(B), mm = (A+B)/2) %>%
    mutate(lnA = log(A), lnB=log(B))
  return(tmp)
}

########################################################
test <- sim.data(25, 3, 0, 1, 0, 0.04, 0.07, 10, 45)   #
view(test)                                             #
########################################################

## Transforms data frames with AR to MOR
## Remember: mm: Mean between measurement methods
## ld: Logarithmic difference
get.mor <- function(clinicals, controls) # Clinical samples and control samples are required to be data frames or data tables
{
  # Transforms AR to MOR
  clinicals <- clinicals %>%
    group_by(sample) %>%
    summarise_at(c("A","B","lnA","lnB","ld","mm"),mean, na.rm=TRUE)
    
  confidence.intervals <- controls %>%
    group_by(sample) %>%
    summarise_at(c("A","B","lnA","lnB","ld","mm"),mean, na.rm=TRUE)
  
  controls <- controls %>%
    mutate_at("sample", as.factor) %>% # Change type of sample column from double to factor (Important when plotting)
    group_by(sample) %>%
    summarise_at(c("A","B","lnA","lnB","ld","mm"),mean, na.rm=TRUE)
    
  return(list(controls = controls, clinicals = clinicals, confidence.intervals = confidence.intervals))
}

## LINEAR ASSUMPTIONS TESTS - DEFAULT - APPROX. 200 TIMES SLOWER THAN "tests.short" ##
get.tests <- function(object,level = 0.95)
{
  o<-object
  pshapiro <- shapiro.test(residuals(o))$p.value
  pbp<- bptest(o)$p.value
  pdw <- dwtest(o)$p.value
  results <- data.table::data.table(N = (pshapiro>(1-level)),
                         H = (pbp>(1-level)),
                         A = (pdw>(1-level))) %>%
    mutate(Total = ifelse(N + H + A == 3, TRUE,FALSE)) %>%
    rename(Normality = "N", Homoscedasticity = "H", Auto.correlation = "A")
  return(results)
}

## LINEAR ASSUMPTIONS CHECKS - "GLVMA" PACKAGE MUST BE LOADED BEFORE USE##
tests.short <- function(glvma.object, level)
{
  ok <- c(glvma.object$GlobalTest$GlobalStat4$pvalue[,],
          glvma.object$GlobalTest$DirectionalStat1$pvalue,
          glvma.object$GlobalTest$DirectionalStat2$pvalue,
          glvma.object$GlobalTest$DirectionalStat3$pvalue[,],
          glvma.object$GlobalTest$DirectionalStat4$pvalue) >= rep(x = 1 - level, 5)
  ifelse(test = sum(ok) == 5, 1, 0)
}

## Automatic checks of commutability (Naive criterion, that is the PI criterion) ##
check.controls <- function(clinicalsAR, clinicals, controls, confidence.intervals)
{
  x <- controls$B # x-values of control material measures. Typically three.
  obj<-lm(data = clinicalsAR, formula = A ~ B) # Linear model
  prediction <- data.table::data.table(predict(object = obj, newdata=list(B=x), level = 0.99, interval = "prediction")) %>%
    dplyr::select(lwr,upr,-fit) %>%
    mutate(mA = confidence.intervals$A, mB = confidence.intervals$B)  %>% ## In mor models we keep the names A,B,ld,mm, etc.
    rowwise() %>%
    summarise(commutable = mA <= upr & mA >= lwr & mB >= min(clinicalsAR$B) & mB <= max(clinicalsAR$B)) # Converts samples that is within prediciton bands 
  return(unname(sum(unlist(prediction))/length(x)))
}

## Check controls of log-log ##
check.controls.loglog <- function(clinicalsAR, clinicals, controls, confidence.intervals)
{
  x <- controls$lnB # log of x of control material measures (MOR)
  new <- get.newdata(clinicals$lnA, clinicals$lnB)
  obj<-lm(data = clinicalsAR, formula = lnA ~ lnB) # log-log model
  prediction <- data.table::data.table(predict(object = obj, newdata=list(lnB=x), level = 0.99, interval = "prediction")) %>%
    dplyr::select(lwr,upr,-fit) %>% # We do not use the fitted coloumn
    mutate(mA = confidence.intervals$lnA, mB = confidence.intervals$lnB) %>%
    rowwise() %>%
    summarise(commutable = mA <= upr & mA >= lwr & mB >= min(clinicalsAR$lnB) & mB <= max(clinicalsAR$lnB))
  prediction<-sum(unlist(prediction))/length(x) # Return the proportion of control material samples being commutable with respect to the naive PI criterion
  return(unname(prediction))
}

## ##
check.controls.ba <- function(clinicalsAR,clinicals,controls,confidence.intervals)
{
  x <- controls$mm # x-values of control material measures (MOR)
  new <- get.newdata(clinicals$mm,clinicals$mm)
  obj<-lm(data=clinicalsAR,formula= ld ~ poly(mm,4)) # BA-model fourth degree
  prediction <- data.table::data.table(predict(object = obj, newdata=list(mm=x), level = 0.99, interval = "prediction")) %>%
    dplyr::select(lwr,upr,-fit) %>%
    mutate(mA = confidence.intervals$ld, mB=confidence.intervals$mm)  %>%
    rowwise() %>%
    summarise(commutable.low = mA <= upr & mA >= lwr & mB >= min(clinicalsAR$mm) & mB <= max(clinicalsAR$mm))
  prediction<-sum(unlist(prediction))/length(x)
  return(unname(prediction))
}

check.controls.dr <- function(clinicalsAR, clinicals, controls, confidence.intervals)
{
  x <- controls$B # x-values of control material measures (MOR)
  new <- get.newdata(clinicals$A,clinicals$B)
  prediction <- deming.predictInterval(clinicalsAR$A,clinicalsAR$B, newdata = x) %>%
    mutate(mA = confidence.intervals$A, mB = confidence.intervals$B) %>% rowwise() %>%
    summarise(commutable = mA <= upr & mA >= lwr & mB >= min(clinicalsAR$B) & mB <= max(clinicalsAR$B))
  return(unname(sum(unlist(prediction))/length(x)))
}

####################################
## Testing #########################
####################################
P <- sim.data(25,3,0,1,0,0.05,0.05,20,60)
C <- sim.data(3,3,0,1,0,0.05,0.05,20,60)

mor <- get.mor(P,C)
morP<-mor$clinicals
morC<-mor$controls
morCI<-mor$confidence.intervals
new <- get.newdata(P$lnA,P$lnB)
obj <- lm(data = P, formula = lnA ~ lnB)
pred <- data.frame(new = new, predict(object = obj, newdata = list(lnB=new), level = 0.99, interval = "prediction"))

check.controls.loglog(clinicalsAR = P, clinicals = morP, controls = morC, confidence.intervals = morCI)
get.tests(obj)$Total

ggplot() + 
  geom_ribbon(data = pred, aes(x = new, ymin = lwr, ymax = upr), alpha = 0.3, fill = "green", color = "black", size = 1) +
  geom_line(data = pred, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_point(data = morP, aes(x = lnB, y = lnA), color = "blue", size = 2) +
  geom_point(data = morC, aes(x = lnB, y = lnA, shape = sample), color = "red", size = 3) +
  xlab("Average of methods") + ylab("Logaritmic difference") + labs(title = "Logaritmic difference vs. average of methods")

###########################################################################
###########################################################################

## DEMING - SIMULATION ##
a<-0 # Non-linearity coefficient
b<-1 # Slope coefficient
c<-0 # Intercept coefficient
r<-3 # Number of replications on samples
n<-25 # Number of patient samples
N<-500 # Number of simulated data-sets.
s<-500 # Number of simulation studies run

## Vectors to be filled in for loop ##
commutability_results <- c()

## Beginning simulation process - Might take some time ##
for (q in 0:s)
{
  ## Simulating the data-sets ##
  simP<-(replicate(N,list(sim.data(n,r,a,b,c,0.3,0.04,70,135))))
  simC<-(replicate(N,list(sim.data(3,r,a,b,c,0.3,0.04*(1+0.01*q),70,135))))
  
  ## Creating a list named "mor" containing all simulated datasets as MOR.
  mor<-list() ## Creating empty list to be filled in for loop ##
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  ## Creating three empty lists to be filled ##
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  
  ## Filling the lists above ##
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  ## Creating two empty vectors to contain acceptance rate of commutability and tests ##
  checks <- c()
  
  ## Filling the empty vectors above ##
  for (i in 1:N)
  {
    checks[i] <- check.controls.dr(clinicalsAR = simP[[i]],clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
  }
  
  ## Saving results as these two variables ##
  commutability_results[q+1] <- (sum(floor(checks)))/N # Rounded down which implies that we are only accepting commutability if all samples are accepted.
}
plot(commutability_results)

different_a_relative <- data.frame(commutability = commutability_results,
                          a = (0:200)/1000) %>%
  gather(key = "Test", value = "Acceptance rate", commutability)

different_a_relative
  
  
plot_different_a_dr <- ggplot(data = different_a_relative,
       aes(x = a, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) + xlab("Non-linearity coefficient relative to patients") +
  labs(title = "Acceptance rates as non-linearity coefficient", subtitle = "relative to patients increases") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)

plot(plot_different_a_dr)

different_b_relative <- data.frame(commutability = commutability_results,
                                   b = (100:200)/100) %>%
  gather(key = "Test", value = "Acceptance rate", commutability)

plot_different_b_dr <- ggplot(data = different_b_relative,
                              aes(x = b, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) + xlab("Slope coefficient relative to patients") +
  labs(title = "Acceptance rates as slope coefficient", subtitle = "relative to patients increases") +
  geom_vline(xintercept = 1, linetype = 2, size = 1)

plot(plot_different_b_dr)

different_c_relative <- data.frame(commutability = commutability_results,
                                   c = (0:1000)/10) %>%
  gather(key = "Test", value = "Acceptance rate", commutability)

plot_different_c_dr <- ggplot(data = different_c_relative,
                              aes(x = c, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2, method = "loess") + xlab("Intercept coefficient relative to patients") +
  labs(title = "Acceptance rates as intercept coefficient", subtitle = "relative to patients increases") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)

plot(plot_different_c_dr)



## Altered acceptance criterion relies on normality of replicates #
p<-sim.data(25,3,0,1,0.001,0.02,0.04,60,110)
c<-sim.data(3,3,0,1,0,0.02,0.04,60,110)

c<-c %>% 
  mutate_at(vars(sample),funs(as.factor)) %>%
  group_by(sample) %>%
  mutate(minA=mean(A) - 2*sd(A),maxA=mean(A) + 2*sd(A),
         minB=mean(B) - 2*sd(B),maxB=mean(B) + 2*sd(B))
  
c <- c %>% mutate(mB=mean(B),mA=mean(A))

predzzzz <-bootstrap.predictInterval(method.A = p$A, method.B = p$B, 5e3, 0.99, 120, 55)

pl <- ggplot() +
  geom_ribbon(data = predzzzz, aes(x=new,ymin=lwr,ymax=upr), fill = "green", alpha = 0.3, color="black", size = 1) +
  geom_point(data = p, aes(x=B,y=A), color = "blue") +
  geom_rect(data = c, aes(xmin=minB,xmax=maxB,ymin=minA,ymax=maxA),fill="red",
            alpha=0.3, size = 1, color = "black") +
  geom_point(data = c, aes(x=mB,y=mA,shape=sample), color = "yellow", size = 3) +
  xlab("Method B") + ylab("Method A") +
  labs(title = "Deming regression assessment procedure",
       subtitle = "With 95% confidence regions")
plot(pl)


## Smoothing splines ##
install.packages(npreg)
require(npreg)


sp<-npreg::sm(data = patients.dim.cob, formula = A ~ B, knots = length(B))
predz<-predict(object = sp, newdata = data.frame(B=patients.dim.cob$B), level=0.99, interval = "prediction")
predz<-data.frame(new=patients.dim.cob$B,predz)
predz

sp2 <- smooth.spline(x=patients.dim.cob$B, y = patients.dim.cob$A)
predz1<-predict(sp2,x=patients.dim.cob$B)

sp2$fit

sp$specs$knots


p1<-ggplot() + 
  geom_ribbon(data = predz, aes(x=new,ymax=upr,ymin=lwr),color="black",fill="green",alpha=0.3) +
  geom_line(data=predz, aes(x=new,y=fit),color="yellow", size=2) +
  geom_line(aes(x=predz1$x, y=predz1$y), size=2,color="violet") +
  geom_point(data = patients.dim.cob, aes(y=A,x=B)) +
  geom_point(data = controls.dim.cob, aes(y=A,x=B, shape=sample), size=4,color = "red") +
  ylab("Method A") + xlab("Method B") + labs(title = "Smoothing splines using npreg::sm(violet) and splines::smooth.spline(yellow)", subtitle = "Prediciton bands created with npreg::sm")
plot(p1)

## Real data analysis ##

## Obtaining and cleaning ##
setwd("~/Masteroppgave filer/Datasett")
df.read.2<-readxl::read_xlsx("clinical_samples_2.xlsx")
clinical.samples.2<-df.read.2 %>%
  slice(1:75,preserve=T) %>%
  slice(-76)%>%
  mutate_all(.funs = as.double)%>%
  mutate_at(vars(sample),.funs=as.factor) %>%
  group_by(sample)
colnames(clinical.samples.2)<-c("sample","replicat","advia","sysmex","cell","abx")
names(clinical.samples.2)

## Thiel-sen regression ##
library(mblm)
ts.model <- mblm(dataframe = clinical.samples.2, formula = sysmex ~ advia)
predeee <- data.frame(new=70:140,predict(object=ts.model, level=0.99, interval = "prediction", newdata = data.frame(advia = 70:140))) 
summary(ts.model)

l <- ggplot() +
  geom_ribbon(data = predeee, aes(x=new,ymin=lwr,ymax=upr), fill="green",
              color = "black", size = 1, alpha = 0.3) +
  geom_line(data = predeee, aes(x=new,y=fit), color="violet", size = 2) +
  geom_point(data = clinical.samples.2, aes(x=advia,y=sysmex)) +
  ylab(names(df.read.2)[4]) + xlab(names(df.read.2)[3]) +
  labs(title = paste(names(df.read.2)[4]," vs. ",names(df.read.2)[3]))


plot(l)

ts.model.2 <- mblm(dataframe = clinical.samples.2, formula = sysmex ~ abx)
predeee <- data.frame(new=70:140,predict(object=ts.model.2, level=0.99, interval = "prediction", newdata = data.frame(abx = 70:140))) 

k <- ggplot() +
  geom_ribbon(data = predeee, aes(x=new,ymin=lwr,ymax=upr), fill="green",
              color = "black", size = 1, alpha = 0.3) +
  geom_line(data = predeee, aes(x=new,y=fit), color="violet", size = 2) +
  geom_point(data = clinical.samples.2, aes(x=abx,y=sysmex)) +
  ylab(names(df.read.2)[4]) + xlab(names(df.read.2)[6]) +
  labs(title = paste(names(df.read.2)[4]," vs. ",names(df.read.2)[6]))


plot(k)

## Smoothing splines models ##
model.sa <- npreg::sm(data = clinical.samples.2, formula = sysmex ~ abx, knots = 75)
pred.sa <- data.table(new=seq(from=77,to=142,by=0.5), predict(object = model.sa, 
                   newdata = data.frame(abx=seq(from=77,to=142,by=0.5)), 
                   level = 0.99, interval = "prediction"))

model.sad <- npreg::sm(data = clinical.samples.2, formula = sysmex ~ advia, knots = 10)
pred.sad <- data.table(new=seq(from=75.3,to=135.3,by=0.5), predict(object = model.sad, 
                                                              newdata = data.frame(advia=seq(from=75.3,to=135.3,by=0.5)), 
                                                              level = 0.99, interval = "prediction"))

model.sc <- npreg::sm(data = clinical.samples.2, formula = sysmex ~ cell,knots = 75)
pred.sc <- data.table(new=seq(from=75.7,to=141.0,by=0.5), predict(object = model.sc, 
            newdata = data.frame(cell=seq(from=75.7,to=141.0,by=0.5)), 
            level = 0.99, interval = "prediction"))

names(df.read.2)
range(clinical.samples.2$cell)

plot.sa <- ggplot() + 
  geom_ribbon(data = pred.sa, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "green", size = 1, alpha = 0.3) +
  geom_line(data = pred.sa, aes(x = new, y = fit), color = "violet" , size = 2) +
  geom_point(data = clinical.samples.2, aes(y = sysmex, x = abx)) +
  ylab(names(df.read.2)[4]) + xlab(names(df.read.2)[6]) + labs(title = "Smoothing splines with clinical samples",
                                                               subtitle = paste(names(df.read.2)[4],"vs.",names(df.read.2)[6]))
plot.sad <- ggplot() + 
  geom_ribbon(data = pred.sad, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "green", size = 1, alpha = 0.3) +
  geom_line(data = pred.sad, aes(x = new, y = fit), color = "violet" , size = 2) +
  geom_point(data = clinical.samples.2, aes(y = sysmex, x = advia)) +
  ylab(names(df.read.2)[4]) + xlab(names(df.read.2)[3]) + labs(title = "Smoothing splines with clinical samples",
                                                               subtitle = paste(names(df.read.2)[4],"vs.",names(df.read.2)[3]))
plot.sc <- ggplot() + 
  geom_ribbon(data = pred.sc, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "green", size = 1, alpha = 0.3) +
  geom_line(data = pred.sc, aes(x = new, y = fit), color = "violet" , size = 2) +
  geom_point(data = clinical.samples.2, aes(y = sysmex, x = cell)) +
  ylab(names(df.read.2)[4]) + xlab(names(df.read.2)[5]) + labs(title = "Smoothing splines with clinical samples",
                                                               subtitle = paste(names(df.read.2)[4],"vs.",names(df.read.2)[5]))
plot(plot.sc)
plot(plot.sa)
plot(plot.sad)

## Using log-log transformation ##
df.read.2log<-readxl::read_xlsx("clinical_samples_2.xlsx")
clinical.samples.2.log <- data.table(df.read.2log) %>%
  slice(-1:-75)%>%
  slice(6:80) %>%
  mutate_all(.funs = as.double)%>%
  mutate_at(vars(sample),.funs=as.factor) %>%
  group_by(sample)
colnames(clinical.samples.2.log)<-c("sample","replicat","advia","sysmex","cell","abx")

model.scl <- npreg::sm(data = clinical.samples.2.log, formula = sysmex ~ cell, knots = 75)
pred.scl <- data.table(new=seq(from=1.879096,to=2.149219,by=0.001), predict(object = model.scl, 
                                                                  newdata = data.frame(cell=seq(from=1.879096,to=2.149219,by=0.001)), 
                                                                  level = 0.99, interval = "prediction"))
pred.scl


plot.scl <- ggplot() + 
  geom_ribbon(data = pred.scl, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "green", size = 1, alpha = 0.3) +
  geom_line(data = pred.scl, aes(x = new, y = fit), color = "violet" , size = 2) +
  geom_point(data = clinical.samples.2.log, aes(y = sysmex, x = cell)) +
  ylab(names(df.read.2)[4]) + xlab(names(df.read.2)[5]) + labs(title = "Smoothing splines with clinical samples",
                                                               subtitle = paste(names(df.read.2)[4],"vs.",names(df.read.2)[5]))

plot(plot.scl)

sim.smooth.data<-sim.data(25,3,0.002,1.2,7,0.04,0.04,70,140)
sim.smooth.controls <- sim.data(3,3,0.002,1.2,8,0.04,0.04,70,140)
sim.smooth.model<- npreg::sm(data = sim.smooth.data, formula = A ~ B, knots = 75)
sim.smooth.predict<-data.frame(new=min(sim.smooth.data$B):max(sim.smooth.data$B),predict(object=sim.smooth.model,
                                                                               newdata = data.frame(B = min(sim.smooth.data$B):max(sim.smooth.data$B)),
                                                                               level = 0.99,
                                                                               interval = "prediction"))


sim.smooth.controls <- sim.smooth.controls %>%
  mutate_at(vars(sample), funs(as.factor))

sim.plot <- ggplot() +
  geom_ribbon(data = sim.smooth.predict, aes(x = new, ymin = lwr, ymax = upr),
              fill = "green", size = 1, color = "black", alpha = 0.3) +
  geom_line(data = sim.smooth.predict, aes(x = new, y = fit), color = "violet", size = 2) +
  geom_point(data = sim.smooth.data, aes(x = B, y = A), color = "blue") +
  geom_point(data = sim.smooth.controls, aes(x = B, y = A, shape = sample), color = "red", size = 3) +
  ylab("Method A") + xlab("Method B") +
  labs(title = "Smoothing splines commutability assessment", subtitle = "Non-linear relationship")

plot(sim.plot)
