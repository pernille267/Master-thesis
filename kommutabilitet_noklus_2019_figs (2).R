if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "splines", "brms", "readxl", "mgcv", 
               "merTools", "gridExtra", "grid", "ggplot2", "microbenchmark")

############# Loading data ###################
setwd("C:/Users/Eier/Downloads")
patients_crp_org <- read_excel('crp_pasient.xlsx')
controls_crp_org <- read_excel('crp_control.xlsx')

################### Data: Architect vs Advia ###########################
patients.arc.adv <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Architect", "Advia") %>%
  na.omit %>%
  rename(A = "Architect", B = "Advia") %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

controls.arc.adv <- dplyr::select(controls_crp_org, sample, replicat, "Architect", "Advia") %>% 
  na.omit %>%
  rename(A = "Architect", B = "Advia") %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

######### Data: Dimension vs Cobas #####################################
patients.dim.cob <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Dimension", "Cobas") %>%
  na.omit %>%
  rename(A = "Dimension", B = "Cobas") %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

controls.dim.cob <- dplyr::select(controls_crp_org, sample, replicat, "Dimension", "Cobas") %>% 
  na.omit %>%
  rename(A = "Dimension", B = "Cobas") %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

######### Data: Advia vs. Dimension ############################################
patients.adv.dim <- patients_crp_org %>%
  dplyr::select(sample, replicat, "Advia", "Dimension") %>%
  na.omit %>%
  rename(A = "Advia", B = "Dimension") %>%
  group_by(sample) %>%
  mutate(mA = mean(A), mB = mean(B)) %>%
  group_by(sample, replicat) %>%
  mutate(mean.bias = (A + B)/2, log.diff = log(A) - log(B))

controls.adv.dim <- dplyr::select(controls_crp_org, sample, replicat, "Advia", "Dimension") %>% 
  na.omit %>%
  rename(A = "Advia", B = "Dimension") %>%
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
min.meas <- min(c(patients.adv.dim$A,patients.adv.dim$B))
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
  geom_point(data = controls.arc.adv, aes(x = B, y = A), color = "red") +
  labs(title = "Scatter plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Architect") +
  xlab("Advia")

# Dimension vs Cobas
plot2b <- ggplot() +
  geom_ribbon(data = pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = intercept.dim.cob, slope = slope.dim.cob, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = B, y = A), color = "red") +
  labs(title = "Scatter plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Dimension") +
  xlab("Cobas")

# Advia vs dimension
plot2c <- ggplot() +
  geom_ribbon(data = pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = intercept.adv.dim, slope = slope.adv.dim, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = B, y = A), color = "red") +
  labs(title = "Standard scatter plot", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Advia") +
  xlab("Dimension")

plot(plot2a)
plot(plot2b)
plot(plot2c)

plot2d <- ggplot() +
  geom_ribbon(data = mean.pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = mean.intercept.arc.adv, slope = mean.slope.arc.adv, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = mB, y = mA), color = "red") +
  labs(title = "Scatter plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Architect") +
  xlab("Advia")

# Dimension vs Cobas
plot2e <- ggplot() +
  geom_ribbon(data = mean.pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = mean.intercept.dim.cob, slope = mean.slope.dim.cob, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = mB, y = mA), color = "red") +
  labs(title = "Scatter plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Dimension") +
  xlab("Cobas")

# Advia vs dimension
plot2f <- ggplot() +
  geom_ribbon(data = mean.pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_abline(intercept = mean.intercept.adv.dim, slope = mean.slope.adv.dim, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = mB, y = mA), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = mB, y = mA), color = "red") +
  labs(title = "Standard scatter plot", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("Advia") +
  xlab("Dimension")

plot(plot2d); plot(plot2e); plot(plot2f) 

#4##### Log-log plots ##########################################################

# Architect vs Advia
plot3a <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "Log-log plot") +
  ylab("log(Architect)") +
  xlab("log(Advia)")

# Dimension vs Cobas
plot3b <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "Log-log plot") +
  ylab("log(Dimension)") +
  xlab("log(Cobas)")

# Advia vs Dimension
plot3c <- ggplot() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "Log-log plot") +
  ylab("log(Advia)") +
  xlab("log(Dimension)")

plot(plot3a)
plot(plot3b)
plot(plot3c)

#5##### Log-log plots with prediction bands ###################################
log.olsr.arc.adv <- lm(data = patients.arc.adv, formula = log(A) ~ log(B))
log.olsr.dim.cob <- lm(data = patients.dim.cob, formula = log(A) ~ log(B))
log.olsr.adv.dim <- lm(data = patients.adv.dim, formula = log(A) ~ log(B))

# Prediction intervals
log.pred.arc.adv <- data.frame(new = log(newdata), predict(object = log.olsr.arc.adv, newdata = list(B = newdata), interval = "prediction", level = 0.99))
log.pred.dim.cob <- data.frame(new = log(newdata), predict(object = log.olsr.dim.cob, newdata = list(B = newdata), interval = "prediction", level = 0.99))
log.pred.adv.dim <- data.frame(new = log(newdata), predict(object = log.olsr.adv.dim, newdata = list(B = newdata), interval = "prediction", level = 0.99))

log.slope.arc.adv <- as.double(olsr.arc.adv$coefficients[2])
log.slope.dim.cob <- as.double(olsr.dim.cob$coefficients[2])
log.slope.adv.dim <- as.double(olsr.adv.dim$coefficients[2])
log.intercept.arc.adv <- as.double(olsr.arc.adv$coefficients[1])
log.intercept.dim.cob <- as.double(olsr.dim.cob$coefficients[1])
log.intercept.adv.dim <- as.double(olsr.adv.dim$coefficients[1])

plot4a <- ggplot() +
  geom_ribbon(data = log.pred.arc.adv, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_line(data = log.pred.arc.adv, aes(x = new, y = fit), color = "black") +
  geom_point(data = patients.arc.adv, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.arc.adv, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "log-log plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("log(Architect)") +
  xlab("log(Advia)")

# Dimension vs Cobas
plot4b <- ggplot() +
  geom_ribbon(data = log.pred.dim.cob, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_line(data = log.pred.dim.cob, aes(x = new, y = fit), color = "black") +
  geom_point(data = patients.dim.cob, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.dim.cob, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "log-log plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("log(Dimension)") +
  xlab("log(Cobas)")

# Advia vs dimension
plot4c <- ggplot() +
  geom_ribbon(data = log.pred.adv.dim, aes(x = new, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black") +
  geom_line(data = log.pred.adv.dim, aes(x = new, y = fit), color = "black") +
  geom_point(data = patients.adv.dim, aes(x = log(B), y = log(A)), color  = "blue") +
  geom_point(data = controls.adv.dim, aes(x = log(B), y = log(A)), color = "red") +
  labs(title = "log-log plot with prediction bands", subtitle = "blue: Clinical samples, red: Control material samples") +
  ylab("log(Advia)") +
  xlab("log(Dimension)")

# Plots
plot(plot4a)
plot(plot4b)
plot(plot4c)

# Comparison between ordinary plots and log-log plots
grid.arrange(plot2a,plot4a, nrow = 1)
grid.arrange(plot2b,plot4b, nrow = 1)
grid.arrange(plot2c,plot4c, nrow = 1)

#6##### Residual plots ########################

plot5a <- ggplot(mapping = aes(x = fitted.values(olsr.arc.adv), y = residuals(olsr.arc.adv)), color = "black") +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  labs(title = "Ordinary plot: Architect vs. Advia") +
  ylab("Residuals") +
  xlab("Fitted values")

plot5b <- ggplot(mapping = aes(x = fitted.values(olsr.dim.cob), y = residuals(olsr.dim.cob)), color = "black") +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  labs(title = "Ordinary plot: Dimension vs. Cobas") +
  ylab("Residuals") +
  xlab("Fitted values")

plot5c <- ggplot(mapping = aes(x = fitted.values(olsr.adv.dim), y = residuals(olsr.adv.dim)), color = "black") +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  labs(title = "Ordinary plot: Advia vs. Dimension") +
  ylab("Residuals") +
  xlab("Fitted values")

plot5d <- ggplot(mapping = aes(x = fitted.values(log.olsr.arc.adv), y = residuals(log.olsr.arc.adv)), color = "black") +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  labs(title = "log-log plot: Architect vs. Advia") +
  ylab("Residuals") +
  xlab("Fitted values")

plot5e <- ggplot(mapping = aes(x = fitted.values(log.olsr.dim.cob), y = residuals(log.olsr.dim.cob)), color = "black") +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  labs(title = "log-log plot: Dimension vs. Cobas") +
  ylab("Residuals") +
  xlab("Fitted values")

plot5f <- ggplot(mapping = aes(x = fitted.values(log.olsr.adv.dim), y = residuals(log.olsr.adv.dim)), color = "black") +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  labs(title = "log-log plot: Advia vs. Dimension") +
  ylab("Residuals") +
  xlab("Fitted values")

plot(plot5a); plot(plot5b); plot(plot5c)
plot(plot5d); plot(plot5e); plot(plot5f)
grid.arrange(plot5a, plot5d, nrow = 1)
grid.arrange(plot5b, plot5e, nrow = 1)
grid.arrange(plot5c, plot5f, nrow = 1)

#7##### Methods of evaluation 2 :  Bland-Altman and polynomial regression (Needs work)
#8###### Deming regression function ############################################

errors.in.variables.lm <- function(method.A = NULL, method.B = NULL, lambda = "u", replicates = 3, omit = FALSE)
{
  ####### Basics ###############################################
  r <- replicates
  n <- length(method.A)/r
  gA <- method.A
  gB <- method.B
  replicat <- sort(rep(1:r, ceiling(n)))
  sample <- rep(1:(ceiling(n)), r)
  ifelse(omit == FALSE, patients.A.B <- data.frame(sample = sample, replicat = replicat, gA = gA, gB = gB), patients.A.B <- data.frame(sample = sample[-omit], replicat = replicat[-omit], gA = gA, gB = gB))
  ##############
  print(patients.A.B)
  anova.A <- lmer(data = patients.A.B, REML = FALSE, formula = gA ~ (1 | sample))
  anova.B <- lmer(data = patients.A.B, REML = FALSE, formula = gB ~ (1 | sample))
  sigma.ee.hat <- as.data.frame(VarCorr(anova.A))[2, 4]
  sigma.uu.hat <- as.data.frame(VarCorr(anova.B))[2, 4]
  lambda.hat <- sigma.ee.hat/sigma.uu.hat
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

deming.lm <- function(method.A, method.B, replicates)
{
  df <- data.table::data.table(sample = rep(1:(length(method.A)/replicates)), replicat = rep(1:replicates, length(method.A)/replicates), A = method.A, B = method.B)
  mean.cov <- data.table::data.table(mA = mean(df$A), mB = mean(df$B), SS.AA = (1/length(df$A))*crossprod(df$A - mean(df$A)), SS.BB = (1/length(df$B))*crossprod(df$B - mean(df$A)),SS.BA = (1/length(df$A))*crossprod(df$A - mean(df$A), df$B - mean(df$B)))
  lambda <- mean.cov$SS.AA / mean.cov$SS.BB
  b1 <- (mean.cov$SS.AA - lambda * mean.cov$SS.BB + sqrt((mean.cov$SS.AA - lambda*mean.cov$SS.BB) ^ 2 + 4 * lambda*mean.cov$SS.BA^2))/(2*mean.cov$SS.BA)
  b0 <- mean.cov$mA - mean.cov$mB * b1
  fitted <- b0 + b1 * df$B
  residuals <- df$A - fitted 
  model <- data.table::data.table(response = df$A, predictor = df$B)
  return(list(residuals = residuals, fitted = fitted, coefficients = data.table::data.table(b0,b1), model.frame = model))
}

#9###### Methods of evaluation 3 : Deming regression prediction function #################################
ccc <- deming.lm(method.A = patients.adv.dim$A, method.B = patients.adv.dim$B, replicates = 3)

d3 <- errors.in.variables.lm(method.A = patients.dim.cob$A, method.B = patients.dim.cob$B)

# Prediction interval for the regression line
errors.in.variables.pi <- function(method.A = NULL, method.B = NULL, newdata = "auto", level = 0.95, replicates = 3, lambda = "u")
{
  if (newdata == "auto") {newdata <- seq(from = min(c(method.A, method.B)), to = max(c(method.A, method.B)), by = 0.01 * max(c(method.A, method.B)))}
  alpha <- 1 - level
  n <- length(method.A)
  fit <- errors.in.variables.lm(method.A = method.A, method.B = method.B)
  lambda <- if (lambda != "u") {lambda} else {fit$lambdas$lambda[1]}
  b0 <- fit$coefficients$deming.coef[1]; b1 <- fit$coefficients$deming.coef[2]
  cov.b0.b1 <- fit$covariance.matrix
  sigma.uu <- fit$sigmas$sigma[2]; sigma.ee <- fit$sigmas$sigma[1]
  t <- qt(1-alpha/2, n-2) # Uppper bound -t will be lower bound from symmetry
  deming.pi.jackknife <- data.frame(newdata) %>%
    mutate(y.pred = b0 + b1 * newdata) %>%
    mutate(pred.se = jackknife.univariate(method.A = method.A, method.B = method.B, fit = newdata)$spot.se.error$std.error) %>%
    mutate(lwr = y.pred - t*pred.se, upr = y.pred + t*pred.se)
  return(deming.pi.jackknife)
}

# Testing prediction interval function
errors.in.variables.pi(method.A = patients.dim.cob$A, method.B = patients.dim.cob$B)
errors.in.variables.pi(method.A = patients.arc.adv$A, method.B = patients.arc.adv$B)
errors.in.variables.pi(method.A = patients.adv.dim$A, method.B = patients.adv.dim$B)

# Time should be improved
system.time(expr = errors.in.variables.pi(method.A = patients.dim.cob$A, method.B = patients.dim.cob$B))

# Is this allowed?
jackknife.univariate <- function(method.A = NULL, method.B = NULL, lambda = "u", estimator = "pred.se", newdata = c(0))
{
  m <- length(method.A) # The number of pairs (A,B)
  n <- length(newdata)
  if (estimator == "pred.se" & n > 0)
  {
    fitt.values <- data.frame(matrix(0, nrow = m, ncol = n)) # data frame to be filled with pred.values
    fitt.values[1,] <- errors.in.variables.lm(method.A = method.A, method.B = method.B, omit = FALSE, fit.vector = fit)$fit.values
    df <- data.frame(A = method.A, B = method.B)
    for (i in 1:(m-1))
    {
      one.out <- data.frame(df[-i,])
      new.fit <- errors.in.variables.lm(method.A = one.out$A, method.B = one.out$B, omit = i, fit.vector = fit)$fit.values
      fitt.values[(i+1),] <- new.fit
    }  
  }
  first.fit <- as.numeric(rep(fitt.values[1,], m))
  first.fit <- matrix(data = first.fit, ncol = n, nrow = m, byrow = TRUE)
  pseudo.var <- m * first.fit - (m-1)*fitt.values
  jack.estimator <- sapply(pseudo.var, FUN = mean)
  jack.matrix <- matrix(data = as.numeric(rep(jack.estimator, m)), ncol = n, nrow = m, byrow = T)
  jack.var <- sapply((1/(m-1))*(pseudo.var - jack.matrix)^2, FUN = sum)
  jack.se <- sapply(jack.var, FUN = sqrt)
  jack.mean.se <- mean(jack.se)
  return(list(jack.std.error.mean = jack.mean.se, spot.se.error = data.frame(row.names = 1:n, std.error = jack.se)))
}  

#################################### Testing the jackknife method ##############
jackknife.univariate(method.A = patients.dim.cob$A, method.B = patients.dim.cob$B, estimator = "pred.se")
jackknife.univariate(method.A = patients.arc.adv$A, method.B = patients.arc.adv$B, estimator = "pred.se")
jackknife.univariate(method.A = patients.adv.dim$A, method.B = patients.adv.dim$B, estimator = "pred.se")
################################################################################


##### PI with bootstrap technique
get.leverage <- function(method.B = NULL)
{
  design.matrix <- as.matrix(data.frame(ones = rep(1, length(method.B)), B = method.B))  
  hat.matrix <- design.matrix%*%solve(t(design.matrix)%*%design.matrix)%*%t(design.matrix)
  leverage <- diag(hat.matrix)
  return(leverage)
}

fit <- deming.lm(method.A = patients.dim.cob$A, method.B = patients.dim.cob$B, replicates = 3)
y.p <- fit$coefficients$b0 + fit$coefficients$b1 * (10:90)
residuals <- fit[[1]]
adjusted.residuals <- (residuals)/(sqrt(1 - get.leverage(patients.dim.cob$B)))
s <- adjusted.residuals - mean(adjusted.residuals)

bootstrap.resample <- function(s, fit, pred.points = 10, replicates = 3)
{
  s <- sample(s,length(s),replace = T)
  y.bs.fit <- unlist(fitted(fit))+s; r <- replicates
  x <- fit$model.frame$predictor
  fit.new <- deming.lm(method.A = y.bs.fit, method.B = x, replicates = r)
  residuals.new <- fit.new$residuals
  adjusted.residuals.new <- (residuals.new)/(sqrt(1 - get.leverage(x)))
  s.bs <- adjusted.residuals.new - mean(adjusted.residuals.new) 
  error.fit <- fit.new$coefficients$b0 - fit.new$coefficients$b0
  error.fit <- error.fit + (fit.new$coefficients$b1 - fit.new$coefficients$b1)*pred.points
  return((unname(error.fit + sample(s.bs, size=1))))
}

replicate(n = 50, expr = bootstrap.resample(s,fit))

draws <- data.table::data.table(matrix(0, nrow = 50, ncol = 81))

for (i in 10:90) {draws[,(i-9)] <- replicate(n = 50, expr = bootstrap.resample(s,fit, pred.points = i))}
for (i in 10:90) {draws[,(i-9)] <- draws[,(i-9)] + y.p[i-9]}
sapply(X = draws + y.p, FUN = quantile, probs = c(0.05,0.95))

hist(draws$V1, breaks = 25)


errors.in.variables.lm(method.A = patients.adv.dim$A,method.B = patients.adv.dim$B)$residuals/sqrt(1-get.leverage(method.A = patients.adv.dim$A, method.B = patients.adv.dim$B)$leverage)


# Prediction intervals for three comparisons
pidr.dim.cob <- errors.in.variables.pi(method.A = patients.dim.cob$A, method.B = patients.dim.cob$B, level = 0.99)
pidr.arc.adv <- errors.in.variables.pi(method.A = patients.arc.adv$A, method.B = patients.arc.adv$B, level = 0.99)
pidr.adv.dim <- errors.in.variables.pi(method.A = patients.adv.dim$A, method.B = patients.adv.dim$B, level = 0.99)

############ Plots of deming with 99% prediction interval ######################
plot6a <- ggplot() +
  geom_ribbon(data = pidr.dim.cob, aes(x = newdata, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = pidr.dim.cob, aes(x = newdata, y = y.pred), size = 0.2, color = "black") +
  geom_point(data = patients.dim.cob, aes(x = B, y = A), color = "blue", alpha = 0.8) +
  geom_point(data = controls.dim.cob, aes(x = B, y = A), color = "red") +
  xlab("Cobas") +
  ylab("Dimension") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

plot6b <- ggplot() +
  geom_ribbon(data = pidr.arc.adv, aes(x = newdata, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = pidr.arc.adv, aes(x = newdata, y = y.pred), size = 0.2, color = "black") +
  geom_point(data = patients.arc.adv, aes(x = B, y = A), color = "blue", alpha = 0.8) +
  geom_point(data = controls.arc.adv, aes(x = B, y = A), color = "red") +
  xlab("Architect") +
  ylab("Advia") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

plot6c <- ggplot() +
  geom_ribbon(data = pidr.adv.dim, aes(x = newdata, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
  geom_line(data = pidr.adv.dim, aes(x = newdata, y = y.pred), size = 0.2, color = "black") +
  geom_point(data = patients.adv.dim, aes(x = B, y = A), color = "blue", alpha = 0.8) +
  geom_point(data = controls.adv.dim, aes(x = B, y = A), color = "red") +
  xlab("Dimension") +
  ylab("Advia") + 
  labs(title = "Deming regression", subtitle = "green region is 99% prediction bands")

# Plots - Deming regression + prediction bands
plot(plot6a)
plot(plot6b)
plot(plot6c)

#10##### Methods of evaluation 4 : Deming regression including variability of control materials
get.confidence.region <- function(method.A = NULL, method.B = NULL, replicates = 3, level = 0.95, normality = TRUE)
{
  alpha <- 1 - level
  r <- replicates
  n <- ceiling(length(method.A) / r)
  if (normality == TRUE)
  {
    t <- qt(1 - alpha / 2, df = r-1)
    df.A.B <- data.frame(sample = rep(1:n, r), replicat = sort(rep(1:r,n)), A = method.A, B = method.B) %>%
      group_by(sample) %>%
      mutate(mA = mean(A), mB = mean(B)) %>%
      mutate(sd.A = sd(A), sd.B = sd(B))
    print(df.A.B)
    
    anova.A <- lmer(data = df.A.B, REML = FALSE, formula = A ~ (1 | sample))
    anova.B <- lmer(data = df.A.B, REML = FALSE, formula = B ~ (1 | sample))
    sigma.ee.hat <- as.data.frame(VarCorr(anova.A))[2, 4]
    sigma.uu.hat <- as.data.frame(VarCorr(anova.B))[2, 4]
    df.A.B <- df.A.B %>% 
      mutate(lwrA = mA - t * (sd.A/sqrt(r)), uprA = mA + t * (sd.A/sqrt(r))) %>%
      mutate(lwrB = mB - t * (sd.B/sqrt(r)), uprB = mB + t * (sd.B/sqrt(r)))
    return(df.A.B[,5:12])
    
  }
  else
  {
    return("hey, still do not know how to deal with lack of normality cases")
    df.A.B <- data.frame(sample = rep(1:n, r), replicat = rep(1:r,n), A = method.A, B = method.B)
    lambda.hat = sigma.ee.hat/sigma.uu.hat
  }
}

# Testing function
get.confidence.region(method.A = controls.dim.cob$A, method.B = controls.dim.cob$B, normality = T)
get.confidence.region(method.A = controls.arc.adv$A, method.B = controls.arc.adv$B, normality = T)
get.confidence.region(method.A = controls.adv.dim$A, method.B = controls.adv.dim$B, normality = T)

# Commutability plot
get.commutability.plot.1 <- function(method.A.controls = NULL, method.B.controls = NULL, method.A.patients = NULL, method.B.patients = NULL, level.pred = 0.99, level.conf = 0.95, replicates = 3, lambda = "u", method.names = c("method A","method B"))
{
  r <- replicates; n <- ceiling(length(method.A.patients)/r); m <- ceiling(length(method.A.controls)/r); N <- n*r; M <- m*r
  pred <- errors.in.variables.pi(method.A = method.A.patients, method.B = method.B.patients)
  df.A.B.controls <- data.frame(sample = rep(1:m,r), replicat = sort(rep(1:r,m)), A = method.A.controls, B = method.B.controls) %>%
    bind_cols(get.confidence.region(method.A = method.A.controls, method.B = method.B.controls, level = level.conf))
  print(df.A.B.controls)
  df.A.B.patients <- data.frame(sample = rep(1:n,r), replicat = sort(rep(1:r,n)), A = method.A.patients, B = method.B.patients)
  print(df.A.B.patients)
  cp <- ggplot() +
    geom_ribbon(data = pred, aes(x = newdata, ymin = lwr, ymax = upr), fill = "green", color = "black", alpha = 0.3) +
    geom_line(data = pred, aes(x = newdata, y = y.pred), size = 0.2, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_point(data = df.A.B.patients, aes(x = B, y = A), color = "blue", alpha = 0.8) +
    geom_rect(data = df.A.B.controls, aes(xmin = lwrB, xmax = uprB, ymin = lwrA, ymax = uprA), fill = "black", color = "black", alpha = 0.1) +
    geom_point(data = df.A.B.controls, aes(x = mB, y = mA), color = "red", shape = 18, size = 3) +
    xlab(method.names[2]) +
    ylab(method.names[1]) + 
    labs(title = "Commutability plot", subtitle = paste(... = "green region is ",as.character(0.99*100),"%"," prediction bands", " and control material samples consist of ", as.character(level.conf * 100),"%"," confidence intervals", sep = ""))
  plot(cp)
}

get.commutability.plot.1(method.A.controls = controls.adv.dim$A, method.B.controls = controls.adv.dim$B, method.A.patients = patients.adv.dim$A, method.B.patients = patients.adv.dim$B, method.names = c("Advia", "Dimension"))
get.commutability.plot.1(method.A.controls = controls.dim.cob$A, method.B.controls = controls.dim.cob$B, method.A.patients = patients.dim.cob$A, method.B.patients = patients.dim.cob$B, method.names = c("Dimension", "Cobas"))
get.commutability.plot.1(method.A.controls = controls.arc.adv$A, method.B.controls = controls.arc.adv$B, method.A.patients = patients.arc.adv$A, method.B.patients = patients.arc.adv$B, method.names = c("Architect", "Advia"))

####### Simulation function ############################################


