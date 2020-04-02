if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "splines", "brms", "readxl", "mgcv", 
               "merTools", "gridExtra", "grid", "ggplot2", "microbenchmark",
               "car", "lme4", "lmtest", "foreach", "doParallel","segmented")
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


#a######## Data: Architect vs Advia ###########################
patients.arc.adv <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Architect", "Advia") %>%
  na.omit %>%
  rename(A = "Architect", B = "Advia") %>%
  mutate_at(vars(sample,replicat), as.factor) %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

controls.arc.adv <- dplyr::select(controls_crp_org, sample, replicat, "Architect", "Advia") %>% 
  na.omit %>%
  rename(A = "Architect", B = "Advia") %>%
  mutate_at(vars(sample,replicat), as.factor) %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

#b######## Data: Dimension vs Cobas #####################################
patients.dim.cob <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Dimension", "Cobas") %>%
  na.omit %>%
  rename(A = "Dimension", B = "Cobas") %>%
  mutate_at(vars(sample,replicat), as.factor) %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

controls.dim.cob <- dplyr::select(controls_crp_org, sample, replicat, "Dimension", "Cobas") %>% 
  na.omit %>%
  rename(A = "Dimension", B = "Cobas") %>%
  mutate_at(vars(sample,replicat), as.factor) %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

#c######## Data: Advia vs. Dimension ############################################
patients.adv.dim <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Advia", "Dimension") %>%
  na.omit %>%
  rename(A = "Advia", B = "Dimension") %>%
  mutate_at(vars(sample,replicat), as.factor) %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

controls.adv.dim <- dplyr::select(controls_crp_org, sample, replicat, "Advia", "Dimension") %>% 
  na.omit %>%
  rename(A = "Advia", B = "Dimension") %>%
  mutate_at(vars(sample,replicat), as.factor) %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))


#1##### Standard scatter plots regarding two measurement methods ###############
# Architect vs Advia
plot1a <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = B, y = A), color = "red") +
  labs(title = "Standard scatter plot", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Architect") +
  xlab("Advia")

# Dimension vs Cobas
plot1b <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = B, y = A), color = "red") +
  labs(title = "Standard scatter plot", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Dimension") +
  xlab("Cobas")

# Advia vs dimension
plot1c <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = B, y = A), color = "red") +
  labs(title = "Standard scatter plot", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Advia") +
  xlab("Dimension")

plot(plot1a)
plot(plot1b)
plot(plot1c)

#2##### Methods of evaluation 1 : Ordinary least squares regression (OLSR) #########
# Linear model objects (all replicates)
olsr.arc.adv <- lm(data = patients.arc.adv, formula = A ~ B)
olsr.dim.cob <- lm(data = patients.dim.cob, formula = A ~ B)
olsr.adv.dim <- lm(data = patients.adv.dim, formula = A ~ B)

# Linear model objects (means of replicates)
mean.olsr.arc.adv <- lm(data = patients.arc.adv, formula = mA ~ mB)
mean.olsr.dim.cob <- lm(data = patients.dim.cob, formula = mA ~ mB)
mean.olsr.adv.dim <- lm(data = patients.adv.dim, formula = mA ~ mB)

# New data
min.meas <- c(min(c(patients.adv.dim$A,patients.adv.dim$B, 3:5)))
max.meas <- max(c(patients.adv.dim$A,patients.adv.dim$B))
newdata <- seq(from = min.meas, to = max.meas, by = 0.5)

# Prediction intervals - All replicates
pred.arc.adv <- data.frame(new = newdata, predict(object = olsr.arc.adv, newdata = list(B = newdata), interval = "prediction", level = 0.99))
pred.dim.cob <- data.frame(new = newdata, predict(object = olsr.dim.cob, newdata = list(B = newdata), interval = "prediction", level = 0.99))
pred.adv.dim <- data.frame(new = newdata, predict(object = olsr.adv.dim, newdata = list(B = newdata), interval = "prediction", level = 0.99))

# Prediction intervals - Mean of replicates
mean.pred.arc.adv <- data.frame(new = newdata, predict(object = mean.olsr.arc.adv, newdata = list(mB = newdata), interval = "prediction", level = 0.99))
mean.pred.dim.cob <- data.frame(new = newdata, predict(object = mean.olsr.dim.cob, newdata = list(mB = newdata), interval = "prediction", level = 0.99))
mean.pred.adv.dim <- data.frame(new = newdata, predict(object = mean.olsr.adv.dim, newdata = list(mB = newdata), interval = "prediction", level = 0.99))


#3##### Plots with 99% prediction bands ################################
# Intercepts ad slopes
slope.arc.adv <- as.double(olsr.arc.adv$coefficients[2])
slope.dim.cob <- as.double(olsr.dim.cob$coefficients[2])
slope.adv.dim <- as.double(olsr.adv.dim$coefficients[2])
mean.slope.arc.adv <-as.double(mean.olsr.arc.adv$coefficients[2])
mean.slope.dim.cob <- as.double(mean.olsr.dim.cob$coefficients[2])
mean.slope.adv.dim <- as.double(mean.olsr.adv.dim$coefficients[2])
intercept.arc.adv <- as.double(olsr.arc.adv$coefficients[1])
intercept.dim.cob <- as.double(olsr.dim.cob$coefficients[1])
intercept.adv.dim <- as.double(olsr.adv.dim$coefficients[1])
mean.intercept.arc.adv <-as.double(mean.olsr.arc.adv$coefficients[1])
mean.intercept.dim.cob <- as.double(mean.olsr.dim.cob$coefficients[1])
mean.intercept.adv.dim <- as.double(mean.olsr.adv.dim$coefficients[1])

plot2a <- ggplot() +
  geom_ribbon(data = pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = intercept.arc.adv, slope = slope.arc.adv, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  labs(title = "Architect vs. Advia", subtitle = "With 99% prediction bands") +
  ylab("Architect") +
  xlab("Advia")

# Dimension vs Cobas
plot2b <- ggplot() +
  geom_ribbon(data = pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = intercept.dim.cob, slope = slope.dim.cob, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  labs(title = "Dimension vs. Cobas", subtitle = "With 99% prediction bands") +
  ylab("Dimension") +
  xlab("Cobas")

# Advia vs dimension
plot2c <- ggplot() +
  geom_ribbon(data = pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = intercept.adv.dim, slope = slope.adv.dim, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  labs(title = "Advia vs. Dimension", subtitle = "With 99% prediction bands") +
  ylab("Advia") +
  xlab("Dimension")

plot(plot2a)
plot(plot2b)
plot(plot2c)

plot2d <- ggplot() +
  geom_ribbon(data = mean.pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = mean.intercept.arc.adv, slope = mean.slope.arc.adv, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = mB, y = mA, shape = sample), size = 3, color = "red") +
  labs(title = "OLSR with MOR", subtitle = "With 99% prediciton bands") +
  ylab("Architect") +
  xlab("Advia")

# Dimension vs Cobas
plot2e <- ggplot() +
  geom_ribbon(data = mean.pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = mean.intercept.dim.cob, slope = mean.slope.dim.cob, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = mB, y = mA, shape = sample), size = 3, color = "red") +
  labs(title = "OLSR with MOR", subtitle = "With 99% prediciton bands") +
  ylab("Dimension") +
  xlab("Cobas")

# Advia vs dimension
plot2f <- ggplot() +
  geom_ribbon(data = mean.pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = mean.intercept.adv.dim, slope = mean.slope.adv.dim, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = mB, y = mA, shape = sample), size = 3, color = "red") +
  labs(title = "OLSR with MOR", subtitle = "Green ribbon is 99% prediciton bands") +
  ylab("Advia") +
  xlab("Dimension")

plot(plot2d); plot(plot2e); plot(plot2f) 

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

#6##### Log-log plots with ols and prediction bands ############################

log.olsr.arc.adv <- lm(data = patients.arc.adv, formula = log(A) ~ log(B))
log.olsr.dim.cob <- lm(data = patients.dim.cob, formula = log(A) ~ log(B))
log.olsr.adv.dim <- lm(data = patients.adv.dim, formula = log(A) ~ log(B))

# Prediction intervals
log.pred.arc.adv <- data.frame(new = log(newdata), predict(object = log.olsr.arc.adv, newdata = list(B = newdata), interval = "prediction", level = 0.99)) 
log.pred.dim.cob <- data.frame(new = log(newdata), predict(object = log.olsr.dim.cob, newdata = list(B = newdata), interval = "prediction", level = 0.99))
log.pred.adv.dim <- data.frame(new = log(newdata), predict(object = log.olsr.adv.dim, newdata = list(B = newdata), interval = "prediction", level = 0.99))

# Coefficients
log.slope.arc.adv <- as.double(olsr.arc.adv$coefficients[2])
log.slope.dim.cob <- as.double(olsr.dim.cob$coefficients[2])
log.slope.adv.dim <- as.double(olsr.adv.dim$coefficients[2])
log.intercept.arc.adv <- as.double(olsr.arc.adv$coefficients[1])
log.intercept.dim.cob <- as.double(olsr.dim.cob$coefficients[1])
log.intercept.adv.dim <- as.double(olsr.adv.dim$coefficients[1])

# Architect vs. Advia
plot4a <- ggplot() +
  geom_ribbon(data = log.pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_line(data = log.pred.arc.adv, aes(x = new, y = fit), color = "black") +
  geom_point(data = patients.arc.adv, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "ln(Architect) vs. ln(Advia)", subtitle = "With 99% prediction bands") +
  ylab("ln(Architect)") + xlab("ln(Advia)")

# Dimension vs Cobas
plot4b <- ggplot() +
  geom_ribbon(data = log.pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_line(data = log.pred.dim.cob, aes(x = new, y = fit), color = "black") +
  geom_point(data = patients.dim.cob, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "ln(Dimension) vs. ln(Cobas)", subtitle = "With 99% prediction bands") +
  ylab("ln(Dimension)") + xlab("ln(Cobas)")

# Advia vs dimension
plot4c <- ggplot() +
  geom_ribbon(data = log.pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_line(data = log.pred.adv.dim, aes(x = new, y = fit), color = "black") +
  geom_point(data = patients.adv.dim, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "ln(Advia) vs. ln(Dimension)", subtitle = "With 99% prediction bands") +
  ylab("ln(Advia)") + xlab("ln(Dimension)")
# Plots
plot(plot4a); plot(plot4b); plot(plot4c)
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

#7##### Methods of evaluation 2 :  Bland-Altman and polynomial regression (Needs work)
#8##### Methods of evaluation 2 : BA-plot + polynomial regression
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

errors.in.variables.lm <- function(method.A = NULL, method.B = NULL, lambda = "u", replicates = 3, omit = FALSE)
{
  ####### Basics ###############################################
  r <- replicates
  n <- length(method.A)/r
  gA <- method.A
  gB <- method.B
  replicat <- sort(rep(1:r, ceiling(n)))
  sample <- rep(1:(ceiling(n)), r)
  ifelse(omit == FALSE, patients.A.B <- data.table::data.table(sample = sample, replicat = replicat, gA = gA, gB = gB), patients.A.B <- data.frame(sample = sample[-omit], replicat = replicat[-omit], gA = gA, gB = gB))
  ##############
  lambda.hat <- estimate.lambda(gA, gB, r)
  lambda.hat. <- ifelse(lambda != "u", lambda, lambda.hat)
  ############# Means and estimated covariances ##################### 
  mA <- mean(gA) # Mean of first group
  mB <- mean(gB) # Mean of second group
  SS.AA <- (1/(n*r))*crossprod(gA-mA)[,] # Vector product, because it is faster
  SS.BB <- (1/(n*r))*crossprod(gB-mB)[,] # Vector product, because it is faster
  SS.BA <- (1/(n*r))*crossprod(gB-mB,gA-mA)[,] # Vector product, because it is fasteer
  
  # Deming coefficients from estimated lambda are always calculated
  A <- sqrt((SS.AA - lambda.hat*SS.BB) ^ 2 + 4 * lambda.hat*SS.BA^2)
  B <- SS.AA - lambda.hat * SS.BB
  b1 <- (B + A) / (2*SS.BA); b0 <- mA - mB*b1
  # If lambda is known, then we get these estimates:
  A. <- sqrt((SS.AA - lambda.hat.*SS.BB) ^ 2 + 4 * lambda.hat.*SS.BA^2)
  B. <- SS.AA - lambda.hat. * SS.BB
  b1. <- (B. + A.) / (2*SS.BA); b0. <- mA - mB*b1.
  ifelse(lambda == "u", b0. <- b0, b0 <- b0)
  ifelse(lambda == "u", b1. <- b1, b0 <- b0)
  # These variances are estimated by using gillard's approach.
  sigma.ll <- (SS.BA)/(b1) # Variation of true value
  sigma.uu <- SS.BB - sigma.ll # Variation of measurement method B
  sigma.ee <- SS.AA - (b1^2) * sigma.ll # Variation of measurement method B
  sigma.det <- SS.BB * SS.AA - SS.BA ^ 2 
  V <- b1^2 * sigma.uu + sigma.ee
  residuals <- gA - (b0 + b1 * gB) 
  
  variance.b1 <- (sigma.det) / (sigma.ll ^ 2)
  variance.b0 <- (mB^2)*((sigma.det)/(sigma.ll ^ 2)) + V
  covariance.b0.b1 <- - mB * variance.b1
  std.error.b1 <- sqrt(abs(variance.b1))
  std.error.b0 <- sqrt(abs(variance.b0))
  cov.mat.b1.bo <- data.table::data.table(b0 = c(variance.b0, covariance.b0.b1), b1 = c(covariance.b0.b1, variance.b1))
  fitted.values <- b0 + b1 * method.B
  
  return(list(
    coefficients = data.table::data.table(Intercept = b0, x = b1),
    altern.coef = data.table::data.table(Intercept = b0., x = b1.),
    lambdas = data.table::data.table(lambda.estimated = lambda.hat, lambda = lambda.hat.),
    covariance.matrix = cov.mat.b1.bo, 
    standard.errors = data.table::data.table(Intercept = std.error.b0, x = std.error.b1), 
    sigmas = data.table::data.table(variable = c("Method A", "Method B", "Latent", "Determinant"), sigma = c(sigma.ee, sigma.uu, sigma.ll, sigma.det)),
    fitted.values = fitted.values,
    residuals = residuals,
    model.frame = data.table::data.table(response = gA, predictor = gB)))
}

### A slightly faster deming lm function - Only essentials ###
deming.lm <- function(method.A, method.B, replicates)
{
  df <- data.table::data.table(sample = rep(1:(length(method.A)/replicates)), replicat = rep(1:replicates, length(method.A)/replicates), A = method.A, B = method.B)
  mean.cov <- data.table::data.table(mA = mean(df$A), mB = mean(df$B), SS.AA = (1/length(df$A))*crossprod(df$A - mean(df$A)), SS.BB = (1/length(df$B))*crossprod(df$B - mean(df$A)),SS.BA = (1/length(df$A))*crossprod(df$A - mean(df$A), df$B - mean(df$B)))
  lambda <- estimate.lambda(method.A = method.A, method.B = method.B)
  b1 <- (mean.cov$SS.AA - lambda * mean.cov$SS.BB + sqrt((mean.cov$SS.AA - lambda*mean.cov$SS.BB) ^ 2 + 4 * lambda*mean.cov$SS.BA^2))/(2*mean.cov$SS.BA)
  b0 <- mean.cov$mA - mean.cov$mB * b1
  fitted <- b0 + b1 * df$B
  residuals <- df$A - fitted 
  model <- data.table::data.table(response = df$A, predictor = df$B)
  return(list(residuals = residuals, fitted = fitted, coefficients = data.table::data.table(b0,b1), model.frame = model))
}

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


## Standard errors of our comparison studies ##
get.se(patients.arc.adv$A, patients.arc.adv$B, 3)
get.se(patients.dim.cob$A,patients.dim.cob$B, 3)
get.se(patients.adv.dim$A, patients.adv.dim$B, 3)

## Coefficients ##
k1<-data.table::data.table(fitted = deming.lm(patients.arc.adv$A, patients.arc.adv$B, 3)$fitted, residuals = deming.lm(patients.arc.adv$A, patients.arc.adv$B, 3)$residuals)
k2<-data.table::data.table(fitted = deming.lm(patients.dim.cob$A, patients.dim.cob$B, 3)$fitted, residuals = deming.lm(patients.dim.cob$A, patients.dim.cob$B, 3)$residuals)
k3<-data.table::data.table(fitted = deming.lm(patients.adv.dim$A, patients.adv.dim$B, 3)$fitted, residuals = deming.lm(patients.adv.dim$A, patients.adv.dim$B, 3)$residuals)

## Residual plots - DR ##

ggplot() + geom_hline(yintercept = 0, size = 1) +
  geom_point(data = k1, aes(x = fitted, y = residuals)) +
  ylab("Residuals") + xlab("Fitted values") +
  labs(title = "Residual vs. fitted plot", subtitle = "Architect vs. Advia")
ggplot() + geom_hline(yintercept = 0, size = 1) +
  geom_point(data = k2, aes(x = fitted, y = residuals)) +
  ylab("Residuals") + xlab("Fitted values") +
  labs(title = "Residual vs. fitted plot", subtitle = "Dimension vs. Cobas")
ggplot() + geom_hline(yintercept = 0, size = 1) +
  geom_point(data = k3, aes(x = fitted, y = residuals)) +
  ylab("Residuals") + xlab("Fitted values") +
  labs(title = "Residual vs. fitted plot", subtitle = "Advia vs. Dimension")












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
  fit <- deming.lm(method.A = method.A, method.B = method.B, replicates = 3)
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
  s <- get.sample(method.A,method.B,3,4)
  fit <- s$fit; s <- s$s; alpha<- 1-level;N<-resamples
  A<-method.A;B<-method.B; new<-lwr:upr
  df<- data.table::data.table(foreach(i=lwr:upr, .packages = c("dplyr", "tidyverse"), .export = c("bootstrap.pred","get.sample", "deming.lm","estimate.lambda","get.leverage", "bootstrap.resample"), .combine=cbind) %dopar% bootstrap.pred(A,B,N,i))
  df<- data.table::data.table(t(sapply(X=df,FUN=quantile,probs=c(alpha/2, 1 - alpha/2)))) %>%
    mutate(pred = fit$coefficients$b0 + fit$coefficients$b1 * new)  
  return(data.table::data.table(new = new, pred = df[,3], lwr=df[,1] + df[,3],upr=df[,2] + df[,3]))
}

pidr.dim.cob <- bootstrap.predictInterval(patients.dim.cob$A,patients.dim.cob$B, 2500, 0.99, 100, 0)
pidr.arc.adv <- bootstrap.predictInterval(patients.arc.adv$A,patients.arc.adv$B, 2500, 0.99, 100, 0)
pidr.adv.dim <- bootstrap.predictInterval(patients.adv.dim$A,patients.adv.dim$B, 2500, 0.99, 100, 0)


#10#### Plots of deming with 99% prediction interval ######################
plot6a <- ggplot() +
  geom_ribbon(data = pidr.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = pidr.dim.cob, aes(x=new,y=pred)) +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color = "blue", alpha = 0.8) +
  geom_point(data = controls.dim.cob, aes(x = B, y = A), color = "red") +
  xlab("Cobas") +
  ylab("Dimension") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

plot6b <- ggplot() +
  geom_ribbon(data = pidr.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = pidr.arc.adv, aes(x = new, y = pred), size = 0.2, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color = "blue", alpha = 0.8) +
  geom_point(data = controls.arc.adv, aes(x = B, y = A), color = "red") +
  xlab("Advia") +
  ylab("Architect") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

plot6c <- ggplot() +
  geom_ribbon(data = pidr.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = pidr.adv.dim, aes(x = new, y = pred), size = 0.2, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color = "blue", alpha = 0.8) +
  geom_point(data = controls.adv.dim, aes(x = B, y = A), color = "red") +
  xlab("Dimension") +
  ylab("Advia") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

# Plots - Deming regression + prediction bands

plot(plot6a)
plot(plot6b)
plot(plot6c)

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
  xlab("Advia") + ylab("Architect") +
  labs(title = "Regression splines assessment evaluation",
       subtitle = "Architect vs. Advia (dashed lines are knots)")
spline.plot.adv.dim <- ggplot() +
  geom_ribbon(data = spline.adv.dim.pred, aes(x = new, ymin = lwr, ymax = upr), fill = "green", color = "black", size = 1, alpha = 0.3) +
  geom_line(data = spline.adv.dim.pred, aes(x = new, y = fit), color = "gray", alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(30, 60), size = 1, linetype = "dashed") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color = "blue") +
  geom_point(data = controls.adv.dim, aes(x = B, y = A, shape = sample), size = 3, color = "red") +
  xlab("Advia") + ylab("Architect") +
  labs(title = "Regression splines assessment evaluation",
       subtitle = "Architect vs. Advia (dashed lines are knots)")

plot(spline.plot.arc.adv)
plot(spline.plot.dim.cob)
plot(spline.plot.adv.dim)

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
  labs(title = "Residual plot - RS", subtitle = "Architect vs. Advia")
res.plot.RS.adv.dim <- ggplot() +
  geom_hline(yintercept = 0, size = 1) +
  geom_point(aes(x = fitted.values(spline.adv.dim.lm), y = residuals(spline.adv.dim.lm))) +
  xlab("Fitted values") + ylab("Residuals") +
  labs(title = "Residual plot - RS", subtitle = "Architect vs. Advia")

plot(res.plot.RS.arc.adv)
plot(res.plot.RS.dim.cob)
plot(res.plot.RS.adv.dim)

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

set.seed(333)

### We are simulating data ###
# can be used both for Patient Samples and Controls, and we can change param between them to test models
sim.data<- function(pairs, replicates, a, b, c, CVX, CVY, lower.limit, upper.limit)
{
  r <- replicates; n <- pairs
  sample <- as.factor(rep(1:n, each = r))
  replicat <- rep(1:r, each = 1, times = n)
  y.true <- rep(runif(n, lower.limit, upper.limit), each = r)
  tmp <- data.table::data.table(sample, replicat, y.true) %>% 
    mutate(x.true = a * y.true ^ 2 + b * y.true + c) %>%
    rowwise() %>%
    mutate(A = y.true * (1 + rnorm(1, 0, CVY))) %>%
    mutate(B = x.true * (1 + rnorm(1, 0, CVX))) %>%
    mutate(ld = log(A) - log(B), mm = (A+B)/2)
  return(tmp)
}


commutability.plot <- function(clinicals, controls, evaluation = "OLSR", level = 0.99)
{
  ev <- evaluation
  clinicals <- clinicals %>%
    dplyr::select(c("sample", "replicat", A, B)) %>%
    drop_na() %>%
    mutate(ld = log(A)-log(B), mm = (A + B)*0.5) %>%
    mutate(lnA = log(A), lnB = log(B))
  
  
  minst <- min(clinicals$A, clinicals$B)
  størst <- max(clinicals$A, clinicals$B)
  minstba <- min(clinicals$mm); størstba <- max(clinicals$mm)
  
  print(c(minstba, størstba))
  print(names(clinicals))
  
  controls <- controls %>%
    dplyr::select(c("sample", "replicat",A, B))%>%
    mutate(ld = log(A)-log(B), mm = (A + B)*0.5) %>%
    drop_na() %>%
    mutate(lnA = log(A), lnB = log(B))
  
  #print(controls)
  
  obj <- lm(formula = A ~ B , data = clinicals)
  obj1 <- lm(formula = log(A) ~ log(B), data = clinicals)
  obj2 <- lm(formula = -ld ~ poly(mm,4), data = clinicals)
  
  pred.olsr <- data.table::data.table(new = (minst:størst*3)/3, predict(object = obj, level=0.99, interval = "prediction", newdata = list(B = (minst:størst*3)/3)))
  pred.ll <- data.table::data.table(new = (minst:størst*3)/3, predict(object = obj1, level=0.99, interval = "prediction", newdata = list(B = (minst:størst*3)/3))) 
  pred.ba <- data.table::data.table(new = (minstba:størstba*3)/3, predict(object = obj2, level=0.99, interval = "prediction", newdata = list(mm = (minstba:størstba*3)/3)))
  olsr <- ggplot() + 
    geom_ribbon(data = pred.olsr, aes(x = new, ymin = lwr, ymax = upr), alpha = 0.3, fill = "green", color = "black", size = 1) +
    geom_smooth(data = clinicals, aes(x = B, y = A), method = "lm", color = "gray", alpha = 0.5) +
    geom_point(data = clinicals, aes(x = B, y = A), color = "blue") +
    geom_point(data = controls, aes(x = B, y = A, shape = sample), color = "red", cex = 3) +
    xlab("Method B") + ylab("method B") + labs(title = "OLSR assessment method")
  
  ll <- ggplot() + geom_ribbon(data = pred.ll, aes(x = log(new), ymin = lwr, ymax = upr), alpha = 0.3, fill = "green", color = "black", size = 1) +
    geom_smooth(data = clinicals, aes(x = lnB, y = lnA), method = "lm", color = "gray", alpha = 0.5) +
    geom_point(data = clinicals, aes(x = lnB, y = lnA), color = "blue") +
    geom_point(data = controls, aes(x = lnB, y = lnA, shape = sample), color = "red", cex = 3) +
    xlab("ln(B)") + ylab("ln(A)") +
    labs(title = "Log-log assessment method")
  
  #ba <- ggplot() 
  
  ifelse(test = ev=="OLSR", yes = plot(olsr), no = ifelse(test = ev=="LL", yes = plot(ll), no = 1))
  
}

commutability.plot(simP.should.ok1, simC.should.ok1, evaluation = "OLSR")





set.seed(3323)

# Simulate clinical samples - Forced linear
simP.should.ok1 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.11, c = 2.4, CVX = 0.02, CVY = 0.04, lower.limit = 5, upper.limit = 90)
simP.should.ok2 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.97, c = -0.8, CVX = 0.03, CVY = 0.02, lower.limit = 5, upper.limit = 90)
simP.should.ok3 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.98, c = -0.7, CVX = 0.05, CVY = 0.01, lower.limit = 5, upper.limit = 90)
simP.should.ok4 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.90, c = -0.1, CVX = 0.03, CVY = 0.01, lower.limit = 5, upper.limit = 90)
simP.should.ok5 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.99, c = 1.2, CVX = 0.01, CVY = 0.05, lower.limit = 5, upper.limit = 90)
simP.should.ok6 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.02, c = 1.9, CVX = 0.05, CVY = 0.03, lower.limit = 5, upper.limit = 90)
simP.should.ok7 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.98, c = -0.2, CVX = 0.03, CVY = 0.07, lower.limit = 5, upper.limit = 90)
simP.should.ok8 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.12, c = -0.2, CVX = 0.01, CVY = 0.05, lower.limit = 5, upper.limit = 90)
simP.should.ok9 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.09, c = -3.1, CVX = 0.07, CVY = 0.08, lower.limit = 5, upper.limit = 90)
simP.should.ok10 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.04, c = -2.2, CVX = 0.03, CVY = 0.1, lower.limit = 5, upper.limit = 90)

# Simulate control samples - Forced linear
simC.should.ok1 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.11, c = 2.4, CVX = 0.02, CVY = 0.04, lower.limit = 5, upper.limit = 90)
simC.should.ok2 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.97, c = -0.8, CVX = 0.03, CVY = 0.02, lower.limit = 5, upper.limit = 90)
simC.should.ok3 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.98, c = -0.7, CVX = 0.05, CVY = 0.01, lower.limit = 5, upper.limit = 90)
simC.should.ok4 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.90, c = -0.1, CVX = 0.03, CVY = 0.01, lower.limit = 5, upper.limit = 90)
simC.should.ok5 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.99, c = 1.2, CVX = 0.01, CVY = 0.05, lower.limit = 5, upper.limit = 90)
simC.should.ok6 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.02, c = 1.9, CVX = 0.05, CVY = 0.03, lower.limit = 5, upper.limit = 90)
simC.should.ok7 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.98, c = -0.2, CVX = 0.03, CVY = 0.07, lower.limit = 5, upper.limit = 90)
simC.should.ok8 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.12, c = -0.2, CVX = 0.01, CVY = 0.05, lower.limit = 5, upper.limit = 90)
simC.should.ok9 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.09, c = -3.1, CVX = 0.07, CVY = 0.08, lower.limit = 5, upper.limit = 90)
simC.should.ok10 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.04, c = -2.2, CVX = 0.03, CVY = 0.1, lower.limit = 5, upper.limit = 90)

### OLSR ###
lmP.should.ok1 <- lm(data = simP.should.ok1, formula = A ~ B)
lmP.should.ok2 <- lm(data = simP.should.ok2, formula = A ~ B)
lmP.should.ok3 <- lm(data = simP.should.ok3, formula = A ~ B)
lmP.should.ok4 <- lm(data = simP.should.ok4, formula = A ~ B)
lmP.should.ok5 <- lm(data = simP.should.ok5, formula = A ~ B)
lmP.should.ok6 <- lm(data = simP.should.ok6, formula = A ~ B)
lmP.should.ok7 <- lm(data = simP.should.ok7, formula = A ~ B)
lmP.should.ok8 <- lm(data = simP.should.ok8, formula = A ~ B)
lmP.should.ok9 <- lm(data = simP.should.ok9, formula = A ~ B)
lmP.should.ok10 <- lm(data = simP.should.ok10, formula = A ~ B)

### OLSR - prediction ###
predP.should.ok1 <- data.table::data.table(new = 5:110, predict(lmP.should.ok1, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok2 <- data.table::data.table(new = 5:110, predict(lmP.should.ok2, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok3 <- data.table::data.table(new = 5:110, predict(lmP.should.ok3, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok4 <- data.table::data.table(new = 5:110, predict(lmP.should.ok4, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok5 <- data.table::data.table(new = 5:110, predict(lmP.should.ok5, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok6 <- data.table::data.table(new = 5:110, predict(lmP.should.ok6, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok7 <- data.table::data.table(new = 5:110, predict(lmP.should.ok7, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok8 <- data.table::data.table(new = 5:110, predict(lmP.should.ok8, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok9 <- data.table::data.table(new = 5:110, predict(lmP.should.ok9, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
predP.should.ok10 <- data.table::data.table(new = 5:110, predict(lmP.should.ok10, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 

### OLSR - Evaluation plots ###
plotP.should.ok1 <- ggplot() +
  geom_ribbon(data = predP.should.ok1, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = predP.should.ok1, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok1, aes(x = B, y = A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok1, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok2 <- ggplot() +
  geom_ribbon(data = predP.should.ok2, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok2, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok2, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok2, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok3 <- ggplot() +
  geom_ribbon(data = predP.should.ok3, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok3, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok3, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok3, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok4 <- ggplot() +
  geom_ribbon(data = predP.should.ok4, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok4, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok4, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok4, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok5 <- ggplot() +
  geom_ribbon(data = predP.should.ok5, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok5, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok5, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok5, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok6 <- ggplot() +
  geom_ribbon(data = predP.should.ok6, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok6, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok6, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok6, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok7 <- ggplot() +
  geom_ribbon(data = predP.should.ok7, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok7, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok7, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok7, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok8 <- ggplot() +
  geom_ribbon(data = predP.should.ok8, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok8, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok8, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok8, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok9 <- ggplot() +
  geom_ribbon(data = predP.should.ok9, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok9, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok9, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok9, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
plotP.should.ok10 <- ggplot() +
  geom_ribbon(data = predP.should.ok10, aes(x = new, ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3) +
  geom_line(data = predP.should.ok10, aes(x = new, y = fit), size = 1) +
  geom_point(data = simP.should.ok10, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok10, aes(x = B, y = A, shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")


plot(plotP.should.ok1);plot(plotP.should.ok2);plot(plotP.should.ok3)
plot(plotP.should.ok4);plot(plotP.should.ok5);plot(plotP.should.ok6)
plot(plotP.should.ok7);plot(plotP.should.ok8);plot(plotP.should.ok9)
plot(plotP.should.ok10)

### normality tests ###
shapiro.test(residuals(lmP.should.ok1));shapiro.test(residuals(lmP.should.ok2));shapiro.test(residuals(lmP.should.ok3))
shapiro.test(residuals(lmP.should.ok4));shapiro.test(residuals(lmP.should.ok5));shapiro.test(residuals(lmP.should.ok6))
shapiro.test(residuals(lmP.should.ok7));shapiro.test(residuals(lmP.should.ok8));shapiro.test(residuals(lmP.should.ok9));shapiro.test(residuals(lmP.should.ok10))

### Almost every test fail regarding normality ###
bptest(lmP.should.ok1);bptest(lmP.should.ok2);bptest(lmP.should.ok3)
bptest(lmP.should.ok4);bptest(lmP.should.ok5)

### Okay we stop there. Every test fail. ### 

### Even though control materials seem to be commutable here ###
### none of the requirements for linear fitting is met. ###
### We can therefore not use any linear model with this ###
### simulated data to get out something trustworthy ###


### OLSR - log-log transformation ###
log.lmP.should.ok1 <- lm(data = simP.should.ok1, formula = log(A) ~ log(B))
log.lmP.should.ok2 <- lm(data = simP.should.ok2, formula = log(A) ~ log(B))
log.lmP.should.ok3 <- lm(data = simP.should.ok3, formula = log(A) ~ log(B))
log.lmP.should.ok4 <- lm(data = simP.should.ok4, formula = log(A) ~ log(B))
log.lmP.should.ok5 <- lm(data = simP.should.ok5, formula = log(A) ~ log(B))
log.lmP.should.ok6 <- lm(data = simP.should.ok6, formula = log(A) ~ log(B))
log.lmP.should.ok7 <- lm(data = simP.should.ok7, formula = log(A) ~ log(B))
log.lmP.should.ok8 <- lm(data = simP.should.ok8, formula = log(A) ~ log(B))
log.lmP.should.ok9 <- lm(data = simP.should.ok9, formula = log(A) ~ log(B))
log.lmP.should.ok10 <- lm(data = simP.should.ok10, formula = log(A) ~ log(B))

### log-log - prediction ###
log.predP.should.ok1 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok1, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok2 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok2, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok3 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok3, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok4 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok4, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok5 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok5, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok6 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok6, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok7 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok7, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok8 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok8, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok9 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok9, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok10 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok10, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 

### OLSR - Evaluation plots ###
log.plotP.should.ok1 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok1, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok1, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok1, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok1, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok2 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok2, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok2, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok2, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok2, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok3 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok3, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok3, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok3, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok3, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok4 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok4, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok4, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok4, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok4, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok5 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok5, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok5, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok5, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok5, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok6 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok6, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok6, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok6, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok6, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok7 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok7, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok7, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok7, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok7, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok8 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok8, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok8, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok8, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok8, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok9 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok9, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok9, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok9, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok9, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok10 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok10, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok10, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok10, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok10, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")


plot(log.plotP.should.ok1);plot(log.plotP.should.ok2);plot(log.plotP.should.ok3)
plot(log.plotP.should.ok4);plot(log.plotP.should.ok5);plot(log.plotP.should.ok6)
plot(log.plotP.should.ok7);plot(log.plotP.should.ok8);plot(log.plotP.should.ok9);plot(log.plotP.should.ok10)


### Normality tests ###
shapiro.test(residuals(log.lmP.should.ok1))
shapiro.test(residuals(log.lmP.should.ok2))
shapiro.test(residuals(log.lmP.should.ok3))
shapiro.test(residuals(log.lmP.should.ok4))
shapiro.test(residuals(log.lmP.should.ok5))
shapiro.test(residuals(log.lmP.should.ok6))
shapiro.test(residuals(log.lmP.should.ok7))
shapiro.test(residuals(log.lmP.should.ok8))
shapiro.test(residuals(log.lmP.should.ok9))
shapiro.test(residuals(log.lmP.should.ok10))

### Equal variance tests ###
bptest(log.lmP.should.ok1)
bptest(log.lmP.should.ok2)
bptest(log.lmP.should.ok3)
bptest(log.lmP.should.ok4)
bptest(log.lmP.should.ok5)
bptest(log.lmP.should.ok6)
bptest(log.lmP.should.ok7)
bptest(log.lmP.should.ok8)
bptest(log.lmP.should.ok9)
bptest(log.lmP.should.ok10)


### BA

### Bland-Altman ###

ba.lmP.should.ok1 <- lm(data = simP.should.ok1, formula = ld ~ poly(mm, 2))
ba.lmP.should.ok2 <- lm(data = simP.should.ok2, formula = ld ~ poly(mm, 2))
ba.lmP.should.ok3 <- lm(data = simP.should.ok3, formula = ld ~ mm)
ba.lmP.should.ok4 <- lm(data = simP.should.ok4, formula = ld ~ mm)
ba.lmP.should.ok5 <- lm(data = simP.should.ok5, formula = ld ~ poly(mm,2))
ba.lmP.should.ok6 <- lm(data = simP.should.ok6, formula = ld ~ poly(mm,2))
ba.lmP.should.ok7 <- lm(data = simP.should.ok7, formula = ld ~ mm)
ba.lmP.should.ok8 <- lm(data = simP.should.ok8, formula = ld ~ mm)
ba.lmP.should.ok9 <- lm(data = simP.should.ok9, formula = ld ~ poly(mm,3))
ba.lmP.should.ok10 <- lm(data = simP.should.ok10, formula = ld ~ mm)


### Predict ###
ba.predP.should.ok1 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok1, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok2 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok2, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok3 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok3, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok4 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok4, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok5 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok5, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok6 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok6, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok7 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok7, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok8 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok8, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok9 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok9, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))
ba.predP.should.ok10 <- data.table::data.table(new = 0:100, predict(ba.lmP.should.ok10, interval = "prediction", level = 0.99, newdata = list(mm = 0:100)))


### BA-plots ###
ba.plotP.should.ok1 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok1, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok1, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok1, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok1, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok2 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok2, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok2, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok2, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok2, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok3 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok3, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok3, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok3, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok3, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok4 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok4, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok4, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok4, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok4, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok5 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok5, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok5, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok5, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok5, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok6 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok6, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok6, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok6, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok6, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok7 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok7, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok7, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok7, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok7, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok8 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok8, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok8, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok8, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok8, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok9 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok9, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok9, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok9, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok9, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)
ba.plotP.should.ok10 <- ggplot() +
  geom_ribbon(data = ba.predP.should.ok10, aes(x = new, ymax = upr, ymin = lwr), color = "black", fill = "yellow", alpha = 0.5, size = 1) +
  geom_line(data = ba.predP.should.ok10, aes(x = new, y = fit), color = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1) + 
  geom_point(data = simP.should.ok10, aes(x = mm, y=ld), color = "blue", size = 2) +
  geom_point(data = simC.should.ok10, aes(x = mm, y=ld, shape = sample), color = "red", size = 3, shape = 17)


plot(ba.plotP.should.ok1);plot(ba.plotP.should.ok2);plot(ba.plotP.should.ok3)
plot(ba.plotP.should.ok4);plot(ba.plotP.should.ok5);plot(ba.plotP.should.ok6)
plot(ba.plotP.should.ok7);plot(ba.plotP.should.ok8);plot(ba.plotP.should.ok9);plot(ba.plotP.should.ok10)

## Formal tests ##
shapiro.test(residuals(ba.lmP.should.ok1))
shapiro.test(residuals(ba.lmP.should.ok2))
shapiro.test(residuals(ba.lmP.should.ok3))
shapiro.test(residuals(ba.lmP.should.ok4))
shapiro.test(residuals(ba.lmP.should.ok5))
shapiro.test(residuals(ba.lmP.should.ok6))
shapiro.test(residuals(ba.lmP.should.ok7))
shapiro.test(residuals(ba.lmP.should.ok8))
shapiro.test(residuals(ba.lmP.should.ok9))
shapiro.test(residuals(ba.lmP.should.ok10))

bptest(ba.lmP.should.ok1)
bptest(ba.lmP.should.ok2)
bptest(ba.lmP.should.ok3)
bptest(ba.lmP.should.ok4)
bptest(ba.lmP.should.ok5)
bptest(ba.lmP.should.ok6)
bptest(ba.lmP.should.ok7)
bptest(ba.lmP.should.ok8)
bptest(ba.lmP.should.ok9)
bptest(ba.lmP.should.ok10)


### Deming procedure ###
pred.fit.dr.should.ok1 <- bootstrap.predictInterval(method.A = simP.should.ok1$A, method.B = simP.should.ok1$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok1$A,simP.should.ok1$B))), lwr = floor(min(pmin(simP.should.ok1$A,simP.should.ok1$B))), level = 0.99)
pred.fit.dr.should.ok2 <- bootstrap.predictInterval(method.A = simP.should.ok2$A, method.B = simP.should.ok2$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok2$A,simP.should.ok2$B))), lwr = floor(min(pmin(simP.should.ok2$A,simP.should.ok2$B))), level = 0.99)
pred.fit.dr.should.ok3 <- bootstrap.predictInterval(method.A = simP.should.ok3$A, method.B = simP.should.ok3$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok3$A,simP.should.ok3$B))), lwr = floor(min(pmin(simP.should.ok3$A,simP.should.ok3$B))), level = 0.99)
pred.fit.dr.should.ok4 <- bootstrap.predictInterval(method.A = simP.should.ok4$A, method.B = simP.should.ok4$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok4$A,simP.should.ok4$B))), lwr = floor(min(pmin(simP.should.ok4$A,simP.should.ok4$B))), level = 0.99)
pred.fit.dr.should.ok5 <- bootstrap.predictInterval(method.A = simP.should.ok5$A, method.B = simP.should.ok5$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok5$A,simP.should.ok5$B))), lwr = floor(min(pmin(simP.should.ok5$A,simP.should.ok5$B))), level = 0.99)
pred.fit.dr.should.ok6 <- bootstrap.predictInterval(method.A = simP.should.ok6$A, method.B = simP.should.ok6$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok6$A,simP.should.ok6$B))), lwr = floor(min(pmin(simP.should.ok6$A,simP.should.ok6$B))), level = 0.99)
pred.fit.dr.should.ok7 <- bootstrap.predictInterval(method.A = simP.should.ok7$A, method.B = simP.should.ok7$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok7$A,simP.should.ok7$B))), lwr = floor(min(pmin(simP.should.ok7$A,simP.should.ok7$B))), level = 0.99)
pred.fit.dr.should.ok8 <- bootstrap.predictInterval(method.A = simP.should.ok8$A, method.B = simP.should.ok8$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok8$A,simP.should.ok8$B))), lwr = floor(min(pmin(simP.should.ok8$A,simP.should.ok8$B))), level = 0.99)
pred.fit.dr.should.ok9 <- bootstrap.predictInterval(method.A = simP.should.ok9$A, method.B = simP.should.ok9$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok9$A,simP.should.ok9$B))), lwr = floor(min(pmin(simP.should.ok9$A,simP.should.ok9$B))), level = 0.99)
pred.fit.dr.should.ok10 <- bootstrap.predictInterval(method.A = simP.should.ok10$A, method.B = simP.should.ok10$B, resamples = 5000, upr = ceiling(max(pmax(simP.should.ok10$A,simP.should.ok10$B))), lwr = floor(min(pmin(simP.should.ok10$A,simP.should.ok10$B))), level = 0.99)


dr.plotP.should.ok1 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok1, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok1, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok1, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok1, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok2 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok2, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok2, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok2, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok2, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok3 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok3, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok3, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok3, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok3, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok4 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok4, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok4, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok4, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok4, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok5 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok5, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok5, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok5, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok5, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok6 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok6, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok6, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok6, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok6, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok7 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok7, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok7, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok7, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok7, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok8 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok8, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok8, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok8, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok8, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok9 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok9, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok9, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok9, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok9, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")
dr.plotP.should.ok10 <- ggplot() +
  geom_ribbon(data = pred.fit.dr.should.ok10, aes(x=new,ymin=lwr,ymax=upr),color="black",fill="yellow",size=1,alpha=0.5) +
  geom_line(data = pred.fit.dr.should.ok10, aes(x=new,y=pred), size = 1) +
  geom_point(data = simP.should.ok10, aes(x=B,y=A), color = "blue", size = 2) +
  geom_point(data = simC.should.ok10, aes(x=B,y=A, shape = sample), color = "red", size = 4) +
  xlab("Measurement method B") + ylab("Measurement method A") +
  labs(title = "A vs. B", subtitle = "Blue: Clinicals; Red: Controls")

plot(dr.plotP.should.ok1);plot(dr.plotP.should.ok2);plot(dr.plotP.should.ok3)
plot(dr.plotP.should.ok4);plot(dr.plotP.should.ok5);plot(dr.plotP.should.ok6)
plot(dr.plotP.should.ok7);plot(dr.plotP.should.ok8);plot(dr.plotP.should.ok9);plot(dr.plotP.should.ok10)

shapiro.test(residuals(deming.lm(simP.should.ok1$A,simP.should.ok1$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok2$A,simP.should.ok2$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok3$A,simP.should.ok3$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok4$A,simP.should.ok4$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok5$A,simP.should.ok5$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok6$A,simP.should.ok6$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok7$A,simP.should.ok7$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok8$A,simP.should.ok8$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok9$A,simP.should.ok9$B,3)))
shapiro.test(residuals(deming.lm(simP.should.ok10$A,simP.should.ok10$B,3)))


## Segmented regression ##
SR.should.ok1 <- segmented.lm(lmP.should.ok1)
SR.should.ok2 <- segmented.lm(lmP.should.ok2)
SR.should.ok3 <- segmented.lm(lmP.should.ok3)
SR.should.ok4 <- segmented.lm(lmP.should.ok4)
SR.should.ok5 <- segmented.lm(lmP.should.ok5)
SR.should.ok6 <- segmented.lm(lmP.should.ok6)
SR.should.ok7 <- segmented.lm(lmP.should.ok7)
SR.should.ok8 <- segmented.lm(lmP.should.ok8)
SR.should.ok9 <- segmented.lm(lmP.should.ok9)
SR.should.ok10 <- segmented.lm(lmP.should.ok10)

## Predictions ##
pred.SR.should.ok1 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok1, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok2 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok2, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok3 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok3, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok4 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok4, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok5 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok5, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok6 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok6, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok7 <-data.table::data.table(new = 0:max(pmax(simP.should.ok7$A,simP.should.ok7$A)), predict.segmented(SR.should.ok7, newdata = data.frame(B = 0:max(pmax(simP.should.ok7$A,simP.should.ok7$A))), level = 0.99, interval = "prediction")) 
pred.SR.should.ok8 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok8, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok9 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok9, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 
pred.SR.should.ok10 <-data.table::data.table(new = 0:100, predict.segmented(SR.should.ok10, newdata = data.frame(B = 0:100), level = 0.99, interval = "prediction")) 


plot.SR.should.ok1 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok1, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok1, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok1, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok2 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok2, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok2, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok2, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok3 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok3, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok3, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok3, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok4 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok4, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok4, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok4, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok5 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok5, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok5, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok5, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok6 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok6, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok6, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok6, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok7 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok7, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok7, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok7, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok8 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok8, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok8, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok8, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok9 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok9, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok9, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok9, aes(x = B, y = A, shape = sample), color = "red", size = 4)
plot.SR.should.ok10 <- ggplot() +
  geom_ribbon(data = pred.SR.should.ok10, aes(x = new, ymin = lwr,ymax = upr), fill = "yellow", color = "black", size =1, alpha = 0.3) +
  geom_point(data = simP.should.ok10, aes(x = B, y = A), color = "blue") +
  geom_point(data = simC.should.ok10, aes(x = B, y = A, shape = sample), color = "red", size = 4)


plot(plot.SR.should.ok1);plot(plot.SR.should.ok2);plot(plot.SR.should.ok3)
plot(plot.SR.should.ok4);plot(plot.SR.should.ok5);plot(plot.SR.should.ok6)
plot(plot.SR.should.ok7);plot(plot.SR.should.ok8);plot(plot.SR.should.ok9)
plot(plot.SR.should.ok10)


shapiro.test(SR.should.ok1$residuals) #
shapiro.test(SR.should.ok2$residuals) #
shapiro.test(SR.should.ok3$residuals)
shapiro.test(SR.should.ok4$residuals) #
shapiro.test(SR.should.ok5$residuals)
shapiro.test(SR.should.ok6$residuals)
shapiro.test(SR.should.ok7$residuals)
shapiro.test(SR.should.ok8$residuals) #
shapiro.test(SR.should.ok9$residuals)
shapiro.test(SR.should.ok10$residuals)

bptest(SR.should.ok1)
bptest(SR.should.ok2)
bptest(SR.should.ok3)
bptest(SR.should.ok4) #
bptest(SR.should.ok5)
bptest(SR.should.ok6)
bptest(SR.should.ok7)
bptest(SR.should.ok8)
bptest(SR.should.ok9)
bptest(SR.should.ok10)

## Larger intercepts ##
set.seed(421) # 1st
set.seed(422) # 2nd
set.seed(423) # 3rd
set.seed(424) # 4th
set.seed(425) # 5th

## Smaller intercepts ##
set.seed(521) # 1st
set.seed(522) # 2nd
set.seed(523) # 3rd
set.seed(524) # 4th
set.seed(525) # 5th


new <- replicate(10, sample(c(-1,1), 1) * runif(1,0,4))


# Simulate clinical samples - Forced linear
simP.should.ok1 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.11, c = new[1], CVX = 0.02, CVY = 0.04, lower.limit = 10, upper.limit = 90)
simP.should.ok2 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.97, c = new[2], CVX = 0.03, CVY = 0.02, lower.limit = 10, upper.limit = 90)
simP.should.ok3 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.98, c = new[3], CVX = 0.05, CVY = 0.01, lower.limit = 10, upper.limit = 90)
simP.should.ok4 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.90, c = new[4], CVX = 0.03, CVY = 0.01, lower.limit = 10, upper.limit = 90)
simP.should.ok5 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.99, c = new[5], CVX = 0.01, CVY = 0.05, lower.limit = 10, upper.limit = 90)
simP.should.ok6 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.02, c = new[6], CVX = 0.05, CVY = 0.03, lower.limit = 10, upper.limit = 90)
simP.should.ok7 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 0.98, c = new[7], CVX = 0.03, CVY = 0.07, lower.limit = 10, upper.limit = 90)
simP.should.ok8 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.12, c = new[8], CVX = 0.01, CVY = 0.05, lower.limit = 10, upper.limit = 90)
simP.should.ok9 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.09, c = new[9], CVX = 0.07, CVY = 0.08, lower.limit = 10, upper.limit = 90)
simP.should.ok10 <- sim.data(pairs = 25, replicates = 3, a = 0, b = 1.04, c = new[10], CVX = 0.03, CVY = 0.1, lower.limit = 10, upper.limit = 90)

# Simulate control samples - Forced linear
simC.should.ok1 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.11, c = new[1], CVX = 0.02, CVY = 0.04, lower.limit = 10, upper.limit = 90)
simC.should.ok2 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.97, c = new[2], CVX = 0.03, CVY = 0.02, lower.limit = 10, upper.limit = 90)
simC.should.ok3 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.98, c = new[3], CVX = 0.05, CVY = 0.01, lower.limit = 10, upper.limit = 90)
simC.should.ok4 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.90, c = new[4], CVX = 0.03, CVY = 0.01, lower.limit = 10, upper.limit = 90)
simC.should.ok5 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.99, c = new[5], CVX = 0.01, CVY = 0.05, lower.limit = 10, upper.limit = 90)
simC.should.ok6 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.02, c = new[6], CVX = 0.05, CVY = 0.03, lower.limit = 10, upper.limit = 90)
simC.should.ok7 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 0.98, c = new[7], CVX = 0.03, CVY = 0.07, lower.limit = 10, upper.limit = 90)
simC.should.ok8 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.12, c = new[8], CVX = 0.01, CVY = 0.05, lower.limit = 10, upper.limit = 90)
simC.should.ok9 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.09, c = new[9], CVX = 0.07, CVY = 0.08, lower.limit = 10, upper.limit = 90)
simC.should.ok10 <- sim.data(pairs = 3, replicates = 3, a = 0, b = 1.04, c = new[10], CVX = 0.03, CVY = 0.1, lower.limit = 10, upper.limit = 90)


### OLSR - log-log transformation ###
log.lmP.should.ok1 <- lm(data = simP.should.ok1, formula = log(A) ~ log(B))
log.lmP.should.ok2 <- lm(data = simP.should.ok2, formula = log(A) ~ log(B))
log.lmP.should.ok3 <- lm(data = simP.should.ok3, formula = log(A) ~ log(B))
log.lmP.should.ok4 <- lm(data = simP.should.ok4, formula = log(A) ~ log(B))
log.lmP.should.ok5 <- lm(data = simP.should.ok5, formula = log(A) ~ log(B))
log.lmP.should.ok6 <- lm(data = simP.should.ok6, formula = log(A) ~ log(B))
log.lmP.should.ok7 <- lm(data = simP.should.ok7, formula = log(A) ~ log(B))
log.lmP.should.ok8 <- lm(data = simP.should.ok8, formula = log(A) ~ log(B))
log.lmP.should.ok9 <- lm(data = simP.should.ok9, formula = log(A) ~ log(B))
log.lmP.should.ok10 <- lm(data = simP.should.ok10, formula = log(A) ~ log(B))

### log-log - prediction ###
log.predP.should.ok1 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok1, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok2 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok2, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok3 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok3, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok4 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok4, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok5 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok5, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok6 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok6, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok7 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok7, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok8 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok8, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok9 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok9, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 
log.predP.should.ok10 <- data.table::data.table(new = 5:110, predict(log.lmP.should.ok10, newdata = list(B = 5:110), level = 0.99, interval = "prediction")) 

### OLSR - Evaluation plots ###
log.plotP.should.ok1 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok1, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok1, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok1, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok1, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok2 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok2, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok2, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok2, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok2, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok3 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok3, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok3, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok3, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok3, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok4 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok4, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok4, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok4, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok4, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok5 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok5, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok5, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok5, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok5, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok6 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok6, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok6, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok6, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok6, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok7 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok7, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok7, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok7, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok7, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok8 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok8, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok8, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok8, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok8, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok9 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok9, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok9, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok9, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok9, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")
log.plotP.should.ok10 <- ggplot() +
  geom_ribbon(data = log.predP.should.ok10, aes(x = log(new), ymin = lwr, ymax = upr), fill = "yellow", color = "black", alpha = .3, size = 1) +
  geom_line(data = log.predP.should.ok10, aes(x = log(new), y = fit), size = 1) +
  geom_point(data = simP.should.ok10, aes(x = log(B), y = log(A)), color = "blue") +
  geom_point(data = simC.should.ok10, aes(x = log(B), y = log(A), shape = sample), color = "red", size = 4) +
  ylab("Measurement method A") + xlab("Measurement method B") +
  labs(title = "Measurement method A vs. Measurement method B")


plot(log.plotP.should.ok1);plot(log.plotP.should.ok2);plot(log.plotP.should.ok3)
plot(log.plotP.should.ok4);plot(log.plotP.should.ok5);plot(log.plotP.should.ok6)
plot(log.plotP.should.ok7);plot(log.plotP.should.ok8);plot(log.plotP.should.ok9);plot(log.plotP.should.ok10)

### Normality 
shapiro.test(residuals(log.lmP.should.ok1)) #
shapiro.test(residuals(log.lmP.should.ok2)) #
shapiro.test(residuals(log.lmP.should.ok3)) 
shapiro.test(residuals(log.lmP.should.ok4)) #
shapiro.test(residuals(log.lmP.should.ok5)) #
shapiro.test(residuals(log.lmP.should.ok6)) #
shapiro.test(residuals(log.lmP.should.ok7)) #
shapiro.test(residuals(log.lmP.should.ok8)) #
shapiro.test(residuals(log.lmP.should.ok9)) #
shapiro.test(residuals(log.lmP.should.ok10)) 

### Equal variance tests ###
bptest(log.lmP.should.ok1) #
bptest(log.lmP.should.ok2) #
bptest(log.lmP.should.ok3) #
bptest(log.lmP.should.ok4) #
bptest(log.lmP.should.ok5) #
bptest(log.lmP.should.ok6) #
bptest(log.lmP.should.ok7) #
bptest(log.lmP.should.ok8) #
bptest(log.lmP.should.ok9) #
bptest(log.lmP.should.ok10) #

###
dwtest(log.lmP.should.ok1) #
dwtest(log.lmP.should.ok2) #
dwtest(log.lmP.should.ok3) 
dwtest(log.lmP.should.ok4) #
dwtest(log.lmP.should.ok5) #
dwtest(log.lmP.should.ok6) 
dwtest(log.lmP.should.ok7) #
dwtest(log.lmP.should.ok8) #
dwtest(log.lmP.should.ok9) #
dwtest(log.lmP.should.ok10) #


print(new)

