## Log-log regression method
data.alt.intercept1 <- c()
data.alt.intercept2 <- c()
data.alt.slope1<-c()
data.alt.slope2<-c()
data.non.lin1<-c()
data.non.lin2<-c()

## a ##
for(q in 0:100)
{
  a<-0;b<-1;c<-0;r<-3;n<-25
  N<-333
  simP<-(replicate(N,list(sim.data(n,r,a,b,c,0.05,0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a + 0.0001*q,b,c,0.05,0.05,30,90))))
  
  mor<-list()
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  checks <- c()
  tests <- c()
  for (i in 1:N)
  {
    checks[i] <- check.controls.loglog(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- get.tests(lm(data=patients[[i]], formula = log(A)~log(B)))$Total
  }
  
  data.non.lin1[q+1] <- (sum(floor(checks)))/N
  data.non.lin2[q+1] <- (sum(floor(tests)))/N
}

## b ##
for(q in 0:100)
{
  a<-0;b<-1;c<-0;r<-3;n<-25
  N<-333
  simP<-(replicate(N,list(sim.data(n,r,a,b,c,0.05,0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a,b+0.01*q,c,0.05,0.05,30,90))))
  
  mor<-list()
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  checks <- c()
  tests <- c()
  for (i in 1:N)
  {
    checks[i] <- check.controls.loglog(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- get.tests(lm(data=patients[[i]], formula = log(A)~log(B)))$Total
  }
  
  data.alt.slope1[q+1] <- (sum(floor(checks)))/N
  data.alt.slope2[q+1] <- (sum(floor(tests)))/N
}

## c ##
for(q in 0:100)
{
  a<-0;b<-1;c<-0;r<-3;n<-25
  N<-333
  simP<-(replicate(N,list(sim.data(n,r,a,b,c,0.05,0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a,b,c+0.2*q,0.05,0.05,30,90))))
  
  mor<-list()
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  checks <- c()
  tests <- c()
  for (i in 1:N)
  {
    checks[i] <- check.controls.loglog(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- get.tests(lm(data=patients[[i]], formula = log(A)~log(B)))$Total
  }
  
  data.alt.intercept1[q+1] <- (sum(floor(checks)))/N
  data.alt.intercept2[q+1] <- (sum(floor(tests)))/N
}

## INTERCEPT ##
alt.int <- data.frame(commutability = data.alt.intercept1,
                      assumptions = data.alt.intercept2,
                      int.alt = seq(from = 0, to = 20, by = 0.2)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)

data.alt.int1

plot.alt.int <- ggplot(data = alt.int, aes(x = int.alt, y = `Acceptance rate`, color = Test)) +
  geom_smooth(size = 2) +
  geom_line(size = 1, alpha = 0.3) + ylab("Acceptance rate") +
  labs(title = "How much alteration of intercept is required", subtitle = "to reject commutability") +
  xlab("Intercept alteration relative to patients") + geom_vline(xintercept = 0, linetype = 2, size = 1)

plot(plot.alt.int)

## SLOPE ##
alt.slope <- data.frame(commutability = data.alt.slope1,
                        assumptions = data.alt.slope2,
                        slope.alt = seq(from = 1, to = 2, by = 0.01)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)


plot.alt.slope <- ggplot(data = alt.slope, aes(x = slope.alt, y = `Acceptance rate`, color = Test)) +
  geom_smooth(size = 2) +
  geom_line(size = 1, alpha = 0.3) + ylab("Acceptance rate") +
  labs(title = "How much alteration of slope is required", subtitle = "to reject commutability") +
  xlab("Alteration of slope relative to patients") + geom_vline(xintercept = 1, linetype = 2, size = 1)
plot(plot.alt.slope)

## NON-LINEARITY ##
non.lin <- data.frame(commutability = data.non.lin1, assumptions = data.non.lin2,
                      non.lin.coef = 0.0001*(0:100)) %>%
  gather(key = "Test", value = "rate", commutability, assumptions)

plot.non.lin <- ggplot(data = non.lin, aes(x = non.lin.coef, y = rate, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + xlab("Non-linear coefficient relative to patients") +
  geom_smooth(size = 2, fill = "gray") +
  ylab("Acceptance rate") + labs(title = "How much non-linearity is required", subtitle = "to reject commutability") + geom_vline(xintercept = 0, linetype = 2, size = 1)

plot(plot.non.lin)

grid.arrange(plot.non.lin,plot.alt.slope,plot.alt.int,nrow =3)


###### OTHER STUFF #######
data.int.magnitude <- c()
data.slope.magnitude <- c()
data.nlin.magnitude <- c()

## a ##
for(q in -100:100)
{
  a<-0;b<-1;c<-0;r<-3;n<-25
  N<-100
  simP<-(replicate(N,list(sim.data(n,r,a,b,c + 0.5*q, 0.05, 0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a,b,c + 0.5*q, 0.05, 0.05,30,90))))
  
  mor<-list()
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  #checks <- c()
  tests <- c()
  for (i in 1:N)
  {
    #checks[i] <- check.controls.loglog(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- get.tests(lm(data=patients[[i]], formula = log(A)~log(B)))$Total
  }
  
  #data.non.lin1[q+1] <- (sum(floor(checks)))/N
  data.int.magnitude[q+101] <- (sum(floor(tests)))/N
}

## b ##
for(q in 0:200)
{
  a<-0;b<-1;c<-0;r<-3;n<-25
  N<-250
  simP<-(replicate(N,list(sim.data(n,r,a,b + 0.05 * q,c, 0.05, 0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a,b + 0.05 * q,c, 0.05, 0.05,30,90))))
  
  mor<-list()
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  #checks <- c()
  tests <- c()
  for (i in 1:N)
  {
    #checks[i] <- check.controls.loglog(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- get.tests(lm(data=patients[[i]], formula = log(A)~log(B)))$Total
  }
  
  #data.non.lin1[q+1] <- (sum(floor(checks)))/N
  data.slope.magnitude[q+1] <- (sum(floor(tests)))/N
}


## c ##
for(q in 0:200)
{
  a<-0;b<-1;c<-0;r<-3;n<-25
  N<-250
  simP<-(replicate(N,list(sim.data(n,r,a+0.01*q,b,c, 0.05, 0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a+0.01*q,b,c, 0.05, 0.05,30,90))))
  
  mor<-list()
  for (i in 1:N)
  {
    mor[[i]] <- get.mor(simP[[i]], simC[[i]])
  }
  
  controls<-list()
  patients<-list()
  confidence.intervals<-list()
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  #checks <- c()
  tests <- c()
  for (i in 1:N)
  {
    #checks[i] <- check.controls.loglog(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- get.tests(lm(data=patients[[i]], formula = log(A)~log(B)))$Total
  }
  
  #data.non.lin1[q+1] <- (sum(floor(checks)))/N
  data.nlin.magnitude[q+1] <- (sum(floor(tests)))/N
}






## Intercept magnitude ##
int.magnitude <- data.frame(assumptions = data.int.magnitude,
                        magnitude = seq(from = -50, to = 50, by = 0.5)) %>%
  gather(key = "Test", value = "Acceptance rate", assumptions)

model <- lm(data = int.magnitude, formula = data.int.magnitude ~ ns(x = magnitude, knots = c(-50,-40,-30, -20,-10,0)))

plot(model$fitted.values)

confint <- data.frame(predict(object = model, interval = "confidence", level = 0.95))

plot.int.magnitude <- ggplot(data = int.magnitude, aes(x = magnitude, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.4) + labs(title = "Linear model assumptions when the intercept alters") +
  geom_ribbon(fill = "gray", color = "gray", aes(x = magnitude, ymin = confint$lwr, ymax = confint$upr)) +
  geom_line(size = 2, aes(x = magnitude, y = model$fitted.values)) +
  xlab("Intercept - c") + geom_vline(xintercept = 0, linetype = 2, size = 1)

plot(plot.int.magnitude)

## Slope magnitude ##
slope.magnitude <- data.frame(assumptions = data.slope.magnitude,
                              magnitude = seq(from = 1, to = 11, by = 0.05)) %>%
  gather(key = "Test", value = "Acceptance rate", assumptions)

plot.slope.magnitude <- ggplot(data = slope.magnitude, aes(x = magnitude, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + labs(title = "Linear model assumptions when  the slope alters") +
  geom_smooth(size = 2) + xlab("Slope coefficient - b") +
  geom_vline(xintercept = 1, size = 1, linetype = 2)

## Non-linearity magnitude ##
nlin.magnitude <- data.frame(assumptions = data.nlin.magnitude,magnitude = seq(from = 0, to = 2, by = 0.01)) %>%
  gather(key = "Test", value = "Acceptance rate", assumptions)

plot.nlin.magnitude <- ggplot(data = nlin.magnitude, aes(x = magnitude, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + labs(title = "Linear model assumptions when the non-linearity coefficient alters") +
  geom_smooth(size = 2) + xlab("Non-lineary coefficient - a") +
  geom_vline(xintercept = 0, size = 1, linetype = 2)


grid.arrange(plot.nlin.magnitude,plot.slope.magnitude,plot.int.magnitude,nrow=3)




