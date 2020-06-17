olsr_ass_b <- c()
olsr_ass_c <- c()

olsr_com_a <- c()
olsr_com_b <- c()
olsr_com_c <- c()

## OLRS method ##
for (q in 0:60)
{
  a<-0;b<-1;c<-0;r<-3
  N<-250
  simP<-(replicate(N,list(sim.data(25,r,a,b,c,0.05,0.05,30,90))))
  simC<-(replicate(N,list(sim.data(3,r,a,b + 0.01*q,c,0.05,0.05,30,90))))
  
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
  #tests <- c()
  for (i in 1:N)
  {
    checks[i] <- check.controls(clinicalsAR = simP[[i]], clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    #tests[i] <- get.tests(lm(data=patients[[i]], formula = A~B))$Total
  }
  
  olsr_com_b[q+1] <- (sum(floor(checks)))/N # Proportion of data-sets fulfilling commutability.
  #olsr_ass_c[q+1] <- (sum(floor(tests)))/N # Proportion of data-sets fulfilling required model assumptions
}

#plot(olsr_diff_n)
#plot(olsr_com_a)
#par(mfrow = c(1,2))
#plot(olsr_com_b)
#plot(olsr_ass_b)
par(mfrow = c(1,1))
plot(olsr_com_c)



## NON-LINEARITY COEFFICENT ##
df_olsr_com_a <- data.frame(commutability = olsr_com_a, a = seq(from = 0, to = 0.006, by = 0.0001)) %>%
  gather(key = "test", value = "Acceptance rate", commutability)
plot_olsr_com_a <- ggplot(data = df_olsr_com_a, aes(x = a, y = `Acceptance rate`,color = test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) + 
  labs(title = "Acceptance concerning alteration of non-linear",
       subtitle = "coefficent deviation betwen controls and patients") +
  geom_vline(size = 1, linetype = 2, xintercept = 0) + xlab("Non-linearity coefficient relative to patients")

plot(plot_olsr_com_a)

## SLOPE COEFFICIENT ##
df_olsr_com_b <- data.frame(commutability = olsr_com_b, b = seq(from = 1, to = 1.6, by = 0.01)) %>%
  gather(key = "test", value = "Acceptance rate", commutability)

plot_olsr_com_b <- ggplot(data = df_olsr_com_b, aes(x = b, y = `Acceptance rate`,color = test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) +
  labs(title = "Acceptance concerning alteration of slope",
       subtitle = "coefficient deviation between controls and patients") +
  geom_vline(size = 1, linetype = 2, xintercept = 1) +
  xlab("Slope coefficient relative to patients")

plot(plot_olsr_com_b)

## INTERCEPT COEFFICIENT ##
df_olsr_com_c <- data.frame(commutability = olsr_com_c, c = seq(from = 0, to = 20, by = 0.2)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability)

plot_olsr_com_c <- ggplot(data = df_olsr_com_c, aes(x = c, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.5) + geom_smooth(size = 2) +
  labs(title = "Acceptance concerning alteration of intercept",
       subtitle = "coefficient deviation between controls and patients") +
  xlab("Intercept coefficient relative to patients") +
  geom_vline(size = 1, linetype = 2, xintercept = 0)

plot(plot_olsr_com_c)

grid.arrange(plot_olsr_com_a,plot_olsr_com_b,plot_olsr_com_c, ncol = 1, nrow = 3)
