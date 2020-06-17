## Commutability checks ##
#ba_com_different_a <- c()
#ba_com_different_b <- c()
ba_com_different_c <- c()
ba_com_different_cvx <- c()

## Linear assumptions checks ##
#ba_ass_different_a <- c()
#ba_ass_different_b <- c()
ba_ass_different_c <- c()
ba_ass_different_cvx <- c()

## BLAND-ALTMAN - SPECIFICATIONs##
a<-0 # Non-linearity coefficient (Common)
b<-1 # Slope coefficient (Common)
c<-0 # Intercept coefficient (Common)
r<-3 # Number of replications on samples (Common)
n<-25 # Number of patient samples (Clinicals)
N<-100 # Number of simulated data-sets. (Common)
s<-100 # Number of simulation studies run (Common)
cvy <- 0.05 # Variability in y-direction (Common)
cvx <- 0.05 # Variability in x-direction (Common)
l <- 30 # Upper range of data
u <- 90 # Lower range of data

## THE SIMULATION SCHEME - REMEMBER TO CHANGE VECTORS ##
for (q in 0:s)
{
  ## SIMULATING N DATA-SETS: BOTH CONTROL AND PATIENTS ##
  simP<-(replicate(N,list(sim.data(n,r,a,b,c,cvx,cvy,l,u))))
  simC<-(replicate(N,list(sim.data(3,r,a,b,c + 0.2 * q,cvx,cvy,l,u))))
  
  ## EMPTY LIST OF MOR DATA FRAMES ##
  mor<-list()
  
  ## FROM AR TO MOR FOR ALL SIMULATED DATA FRAMES ##
  for (i in 1:N){mor[[i]] <- get.mor(simP[[i]], simC[[i]])}
  
  ## MAKING THREE LISTS TO CONTAIN EACH COMPONENT OF MOR OUTPUT ##
  controls<-list(); patients<-list(); confidence.intervals<-list()
  
  ## Filling the lists above ##
  for (i in 1:N)
  {
    controls[[i]] <- mor[[i]]$controls
    patients[[i]] <- mor[[i]]$clinicals
    confidence.intervals[[i]] <- mor[[i]]$confidence.intervals
  }
  
  ## RESULTS ARE SAVED IN THESE VECTORS ##
  checks <- c(); tests <- c()
  
  ## Filling the empty vectors above ##
  for (i in 1:N)
  {
    checks[i] <- check.controls.ba(clinicalsAR = simP[[i]],clinicals = patients[[i]], controls = controls[[i]], confidence.intervals = confidence.intervals[[i]])
    tests[i] <- tests.short(glvma.object = gvlma.lm(lm(data = simP[[i]], formula = ld ~ poly(mm,4))),level=0.95)
  }
  
  ## 1 - ALL OR NONE ##
  ## UNCOMMENT THE ONES YOU  WANT TO FILL - REMEMBER TO SPECIFY THE INDEX CORRECTLY##
  
  ## 1. Different non-linearity ##
  #ba_com_different_a[q+1] <- (sum(floor(checks)))/N
  #ba_ass_different_a[q+1] <- (sum(floor(tests)))/N
  
  ## 2. Different slope ##
  #ba_com_different_b[q+1] <- (sum(floor(checks)))/N
  #ba_ass_different_b[q+1] <- (sum(floor(tests)))/N
  
  ## 3. Diferent intercept ##
  ba_com_different_c[q+1] <- (sum(floor(checks)))/N
  ba_ass_different_c[q+1] <- (sum(floor(tests)))/N
  
  ## 4. Different variability in x ##
  #ba_com_different_cvx[q+1] <- (sum(floor(checks)))/N
  #ba_ass_different_cvx[q+1] <- (sum(floor(tests)))/N
}


## What happens when n (sample size) increases ##
different_n <- data.frame(commutability = commutability_results,
                          assumptions = assumptions_restuls,
                          sample.size = 25:100) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)

plot_different_n <- ggplot(data = different_n, aes(x = sample.size, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1) + geom_smooth(size = 2) +
  xlab("Sample size (n)") +
  labs(title = "Acceptance rate when n increases")

plot(plot_different_n)

## What happens when the number of replications increases ##
different_r <- data.frame(commutability = commutability_results,
                          assumptions = assumptions_restuls,
                          replicates = 2:20) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)

plot_different_r <- ggplot(data = different_r, aes(x = replicates, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1) + geom_smooth(size = 2) +
  xlab("Number of replications on every sample (r)") +
  labs(title = "Acceptance rate when r increases")

plot(plot_different_r)

## Plots together ##
grid.arrange(plot_different_n,plot_different_r,nrow=2)

## What happens when non-linear coefficient, slope coefficient and intercept coefficient differs from their typical values ##
different_a <- data.frame(commutability = ba_com_different_a,
                          assumptions = ba_ass_different_a,
                          a = seq(from = 0, to = 0.01, by = 0.0001)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)

plot_different_a <- ggplot(data = different_a,
                           aes(x = a, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) + xlab("Non-linearity coefficient relative to patients") +
  labs(title = "Acceptance rates as non-linearity coefficient", subtitle = "relative to patients increases") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)
plot(plot_different_a)

different_b <- data.frame(commutability = ba_com_different_b,
                          assumptions = ba_ass_different_b,
                          b = seq(from=1,to=1.5,by=0.005)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)

plot_different_b <- ggplot(data = different_b,
                           aes(x = b, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) + xlab("Slope coefficient relative to patients") +
  labs(title = "Acceptance rates as slope coefficient", subtitle = "relative to patients increases") +
  geom_vline(xintercept = 1, linetype = 2, size = 1)
plot(plot_different_b)

different_c <- data.frame(commutability = ba_com_different_c,
                          assumptions = ba_ass_different_c,
                          c = seq(from=0,to=20,by=0.2)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)

plot_different_c <- ggplot(data = different_c,
                           aes(x = c, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) + xlab("Slope coefficient relative to patients") +
  labs(title = "Acceptance rates as intercept coefficient", subtitle = "relative to patients increases") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)
plot(plot_different_c)

different_cvx <- data.frame(commutability = commutability_results,
                            assumptions = assumptions_restuls,
                            CVX = seq(from=0,to=0.2,by=0.002)) %>%
  gather(key = "Test", value = "Acceptance rate", commutability, assumptions)


plot_different_cvx <- ggplot(data = different_cvx,
                             aes(x=CVX, y = `Acceptance rate`, color = Test)) +
  geom_line(size = 1, alpha = 0.3) + geom_smooth(size = 2) +
  geom_vline(xintercept = 0.04, linetype = 2) + xlab("Coefficient of variation in x-direction") +
  labs(title = "Effect of ignored variability in x")
plot(plot_different_cvx)

grid.arrange(plot_different_a,plot_different_b,plot_different_c)

