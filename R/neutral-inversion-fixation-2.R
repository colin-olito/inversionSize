# Data storage
k.vals     <-  c()
inv.len    <-  c()
results    <-  c()

# Parameters to vary
Ns  <-  c(10^4, 10^5) #, 10^6) Simulations for 10^6 take a LONG time!!!
Us  <-  c(0.02, 0.05, 0.1)
s.d  <-  0.01

# Constant parameters
inv.len  <- c(1:10)/10 - 0.05 #inversion lengths

# loop over Ns
for(n in 1:length(Ns)) {
  print(paste("N = ", Ns[n]))
  N.0  <-  Ns[n]
  N    <-  N.0
  q.0  <-  1/N

  # loop over selection strength
  for (u in 1:length(Us)) {
    print(paste("U = ", Us[u]))
    U  <-  Us[u]
    
    #distribution of fitness effects
    # Small k
    k     = 2
    print(paste("k = ", k))
    theta = s.d/k
    s.d.i = rgamma(10^7, k, scale = theta) #1 million selection coeff.
    s.H   = 1/mean(1/s.d.i) #harmonic mean deleterious selection coefficient

    #number of simulations per inversion length
    sims        = 100*N
    sim.results = rep(0, 10)

    # loop over inversion lengths
    for(l in 1:10){
      x = l/10 - 0.05 #inversion length
  
      #cumulative fixations
      sum.fix = 0
    
      # loop over replicate simulations
      for(j in 1:sims){
        q = q.0
        t = 0
    
        #simulation runs as long as there is variation and Ns > 1
        while((N*(exp(U*x*(1 + theta*t)^(-k))-1) > 1) & (q*(1 - q) > 0)){
          #expected frequency after selection
          w.avg = q*exp(-U*x*(1 - (1 + theta*t)^(-k))) + (1 - q)*exp(-U*x)
          F.sel = q + q*(1 - q)*(exp(-U*x*(1 - (1 + theta*t)^(-k))) - exp(-U*x))/w.avg
      
          #binomial sampling
          q = rbinom(1, N, F.sel)/(N)
          t = t + 1
        }
    
        #sum fixations
        sum.fix = sum.fix + rbinom(1, 1, q)
        cat('\r', paste("Progress:", round(100*(j/sims)), "% complete"))
      }
  
      sim.results[l] = exp(-U*x/s.H)*sum.fix/sims
      print(l)
    }
  results  <-  c(results, sim.results)
  inv.len  <-  c(inv.len, x)
  k.vals   <-  c(k.vals, rep(k, times=length(sim.results)))


    #intermediate k simulations
    #distribution of fitness effects
    k = 3
    print(paste("k = ", k))
    theta = s.d/k
    s.d.i = rgamma(10^7, k, scale = theta) #1 million selection coeff.
    s.H = 1/mean(1/s.d.i) #harmonic mean deleterious selection coefficient
    q.0 = 1/N
    
    #number of simulations per inversion length
    sim.results = rep(0, 10)
    
    # loop over inversion lengths
    for(l in 1:10){
      x = l/10 - 0.05 #inversion length
      
      #cumulative fixations
      sum.fix = 0
      
      # loop over replicate simulations
      for(j in 1:sims){
        q = q.0
        t = 0
        
        #simulation runs as long as there is variation and Ns > 1
        while((N*(exp(U*x*(1 + theta*t)^(-k))-1) > 1) & (q*(1 - q) > 0)){
          #expected frequency after selection
          w.avg = q*exp(-U*x*(1 - (1 + theta*t)^(-k))) + (1 - q)*exp(-U*x)
          F.sel = q + q*(1 - q)*(exp(-U*x*(1 - (1 + theta*t)^(-k))) - exp(-U*x))/w.avg
          
          #binomial sampling
          q = rbinom(1, N, F.sel)/(N)
          t = t + 1
        }
        
        #sum fixations
        sum.fix = sum.fix + rbinom(1, 1, q)
        cat('\r', paste("Progress:", round(100*(j/sims)), "% complete"))
      }
      
      sim.results[l] = exp(-U*x/s.H)*sum.fix/sims
      print(l)
    }
  results  <-  c(results, sim.results)
  k.vals  <-  c(k.vals,rep(k, times=length(sim.results)))


    #large k simulations
    #distribution of fitness effects
    k = 4
    print(paste("k = ", k))
    theta = s.d/k
    s.d.i = rgamma(10^7, k, scale = theta) #1 million selection coeff.
    s.H = 1/mean(1/s.d.i) #harmonic mean deleterious selection coefficient
    q.0 = 1/N
    
    #number of simulations per inversion length
    sim.results = rep(0, 10)
      
    # loop over inversion lengths
    for(l in 1:10){
      x = l/10 - 0.05 #inversion length
      
      #cumulative fixations
      sum.fix = 0
      
      # loop over replicate simulations
      for(j in 1:sims){
        q = q.0
        t = 0
        
        #simulation runs as long as there is variation and Ns > 1
        while((N*(exp(U*x*(1 + theta*t)^(-k))-1) > 1) & (q*(1 - q) > 0)){
          #expected frequency after selection
          w.avg = q*exp(-U*x*(1 - (1 + theta*t)^(-k))) + (1 - q)*exp(-U*x)
          F.sel = q + q*(1 - q)*(exp(-U*x*(1 - (1 + theta*t)^(-k))) - exp(-U*x))/w.avg
          
          #binomial sampling
          q = rbinom(1, N, F.sel)/(N)
          t = t + 1
        }
        
        #sum fixations
        sum.fix = sum.fix + rbinom(1, 1, q)
        cat('\r', paste("Progress:", round(100*(j/sims)), "% complete"))
      }
      
      sim.results[l] = exp(-U*x/s.H)*sum.fix/sims
      print(l)
    }
  results  <-  c(results, sim.results)
  k.vals  <-  c(k.vals,rep(k, times=length(sim.results)))

    #######################
    #fixed s.d simulations
    print(paste("fixed s.d"))
    q.0 = 1/N
    
    #number of simulations per inversion length
    sim.results = rep(0, 10)

    # loop over inversion lengths    
    for(l in 1:10){
      x = l/10 - 0.05 #inversion length
      
      #cumulative fixations
      sum.fix = 0
      
      # loop over replicate simulations
      for(j in 1:sims){
        q = q.0
        t = 0
        
        #simulation runs as long as there is variation and Ns > 1
        while((N*(exp(U*x*exp(-s.d*t))-1) > 1) & (q*(1 - q) > 0)){
          #expected frequency after selection
          w.avg = q*exp(-U*x*(1 - exp(-s.d*t))) + (1 - q)*exp(-U*x)
          F.sel = q + q*(1 - q)*(exp(-U*x*(1 - exp(-s.d*t))) - exp(-U*x))/w.avg
          
          #binomial sampling
          q = rbinom(1, N, F.sel)/(N)
          t = t + 1
        }
        
        #sum fixations
        sum.fix = sum.fix + rbinom(1, 1, q)
        cat('\r', paste("Progress:", round(100*(j/sims)), "% complete"))
      }
      
      sim.results[l] = exp(-U*x/s.d)*sum.fix/sims
      print(l)
    }
  results  <-  c(results, sim.results)
  k.vals  <-  c(k.vals,rep(NA, times=length(sim.results)))
  
  }

}

Udf   <-  rep(rep(Us, each = 4*length(inv.len)), times=length(Ns))
NsDf  <-  rep(Ns, each = length(Us)*4*length(inv.len))
invL  <-  rep(inv.len, times=length(Us)*4*length(Ns))
sds   <-  rep(s.d, times=length(results))
out.dat  <-  data.frame(cbind(NsDf,Udf,k.vals, sds, invL, results))

colnames (out.dat)  <-  c("N", "sd", "k", "invLen", "Pfix")
write.csv(out.dat, file="neutralPfixData_U0.2.csv", row.names=FALSE)


dat  <-  read.csv(file="./neutralPfixData_U0.05.csv", header=TRUE)
dat$Pfix[dat$N == 1e+04]
dat$Pfix[dat$N == 1e+05]
N1k  <-  dat[dat$N == 1e+04,]
N10k  <-  dat[dat$N == 1e+05,]

pdf("neutralFixDFE_U0.05.pdf", 7, 7)
par(mfrow = c(2, 2))
plot(0:10/10, rep(1/N1k$N[1],11), type = "l", ylim = c(0, 1.25/N1k$N[1]), xlab = "inversion size (x)", ylab = "fixation probability", main = "s = 0.02; N = 10,000", pch = 19)
points(Pfix[N1k$s == 0.02 & N1k$k == 2] ~ 
     invLen[N1k$s == 0.02 & N1k$k == 2], pch = 21, bg = "orange", data=N1k)
points(Pfix[N1k$s == 0.02 & N1k$k == 3] ~ 
     invLen[N1k$s == 0.02 & N1k$k == 3], pch = 21, bg = "skyblue", data=N1k)
points(Pfix[N1k$s == 0.02 & N1k$k == 4] ~ 
     invLen[N1k$s == 0.02 & N1k$k == 4], pch = 21, bg = "green", data=N1k)
points(Pfix[N1k$s == 0.02 & is.na(N1k$k) == TRUE] ~ 
     invLen[N1k$s == 0.02 & is.na(N1k$k) == TRUE], pch = 21, bg = "white", data=N1k)

plot(0:10/10, rep(1/N10k$N[1],11), type = "l", ylim = c(0, 1.25/N10k$N[1]), xlab = "inversion size (x)", ylab = "fixation probability", main = "s = 0.02; N = 100,000", pch = 19)
points(Pfix[N10k$s == 0.02 & N10k$k == 2] ~ 
     invLen[N10k$s == 0.02 & N10k$k == 2], pch = 21, bg = "orange", data=N10k)
points(Pfix[N10k$s == 0.02 & N10k$k == 3] ~ 
     invLen[N10k$s == 0.02 & N10k$k == 3], pch = 21, bg = "skyblue", data=N10k)
points(Pfix[N10k$s == 0.02 & N10k$k == 4] ~ 
     invLen[N10k$s == 0.02 & N10k$k == 4], pch = 21, bg = "green", data=N10k)
points(Pfix[N10k$s == 0.02 & is.na(N10k$k) == TRUE] ~ 
     invLen[N10k$s == 0.02 & is.na(N10k$k) == TRUE], pch = 21, bg = "white", data=N10k)
legend( x=0.05, y=6e-06,
       legend = c("k = 2",
                  "k = 3",
                  "k = 4",
                  expression(paste("Constant ", s[d]))),
       pch   = 21,
       pt.bg = c("orange", "skyblue", "green", "white"),
       bty     =  'n',
       border  =  NA  
        )
plot(0:10/10, rep(1/N1k$N[1],11), type = "l", ylim = c(0, 1.25/N1k$N[1]), xlab = "inversion size (x)", ylab = "fixation probability", main = "s = 0.05; N = 10,000", pch = 19)
points(Pfix[N1k$s == 0.05 & N1k$k == 2] ~ 
     invLen[N1k$s == 0.05 & N1k$k == 2], pch = 21, bg = "orange", data=N1k)
points(Pfix[N1k$s == 0.05 & N1k$k == 3] ~ 
     invLen[N1k$s == 0.05 & N1k$k == 3], pch = 21, bg = "skyblue", data=N1k)
points(Pfix[N1k$s == 0.05 & N1k$k == 4] ~ 
     invLen[N1k$s == 0.05 & N1k$k == 4], pch = 21, bg = "green", data=N1k)
points(Pfix[N1k$s == 0.05 & is.na(N1k$k) == TRUE] ~ 
     invLen[N1k$s == 0.05 & is.na(N1k$k) == TRUE], pch = 21, bg = "white", data=N1k)

plot(0:10/10, rep(1/N10k$N[1],11), type = "l", ylim = c(0, 1.25/N10k$N[1]), xlab = "inversion size (x)", ylab = "fixation probability", main = "s = 0.05; N = 100,000", pch = 19)
points(Pfix[N10k$s == 0.05 & N10k$k == 2] ~ 
     invLen[N10k$s == 0.05 & N10k$k == 2], pch = 21, bg = "orange", data=N10k)
points(Pfix[N10k$s == 0.05 & N10k$k == 3] ~ 
     invLen[N10k$s == 0.05 & N10k$k == 3], pch = 21, bg = "skyblue", data=N10k)
points(Pfix[N10k$s == 0.05 & N10k$k == 4] ~ 
     invLen[N10k$s == 0.05 & N10k$k == 4], pch = 21, bg = "green", data=N10k)
points(Pfix[N10k$s == 0.05 & is.na(N10k$k) == TRUE] ~ 
     invLen[N10k$s == 0.05 & is.na(N10k$k) == TRUE], pch = 21, bg = "white", data=N10k)
dev.off()