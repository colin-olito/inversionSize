U = 0.05
s.d = 0.01
N.0 = 10^3
N = N.0
q.0 = 1/N

#number of simulations per inversion length
sims = 100*N
sim.results = rep(0, 10)

for(k in 1:10){
  x = k/10 - 0.05 #inversion length
  
  #number of generations per simulation (time until inversion is effectively neutral Ns_t = 1)
  gens = round(max(1, (1/s.d)*log(N*U*x)))
  
  #cumulative fixations
  sum.fix = 0
  
  for(j in 1:sims){
    q = q.0
    t = 0
    
    while((t < gens) & (q > 0)){
      #expected frequency after selection
      w.avg = q*exp(-U*x*(1 - exp(-s.d*t))) + (1 - q)*exp(-U*x)
      F.sel = q + q*(1 - q)*(exp(-U*x*(1 - exp(-s.d*t))) - exp(-U*x))/w.avg
      
      #binomial sampling
      q = rbinom(1, N, F.sel)/(N)
      t = t + 1
    }
    
    #sum fixations
    sum.fix = sum.fix + rbinom(1, 1, q)
  }
  
  sim.results[k] = exp(-U*x/s.d)*sum.fix/sims
  print(k)
}

small.N = sim.results

par(mfrow = c(1, 1))
plot(0:10/10, rep(1/N.0,11), type = "l", ylim = c(0, 1.25/N.0), xlab = "inversion size (x)", ylab = "fixation probability", pch = 19)
points((1:10/10 - 0.05), small.N, pch = 21, bg = "WHITE")


#intermediate N simulations
N = N.0*10
q.0 = 1/N

#number of simulations per inversion length
sims = 100*N
sim.results = rep(0, 10)

for(k in 1:10){
  x = k/10 - 0.05 #inversion length
  
  #number of generations per simulation (time until inversion is effectively neutral Ns_t = 1)
  gens = round(max(1, (1/s.d)*log(N*U*x)))
  
  #cumulative fixations
  sum.fix = 0
  
  for(j in 1:sims){
    q = q.0
    t = 0
    
    while((t < gens) & (q > 0)){
      #expected frequency after selection
      w.avg = q*exp(-U*x*(1 - exp(-s.d*t))) + (1 - q)*exp(-U*x)
      F.sel = q + q*(1 - q)*(exp(-U*x*(1 - exp(-s.d*t))) - exp(-U*x))/w.avg
      
      #binomial sampling
      q = rbinom(1, N, F.sel)/(N)
      t = t + 1
    }
    
    #sum fixations
    sum.fix = sum.fix + rbinom(1, 1, q)
  }
  
  sim.results[k] = exp(-U*x/s.d)*sum.fix/sims
  print(k)
}

inter.N = sim.results

points((1:10/10 - 0.05),inter.N*10, pch = 21, bg = "GRAY")


#large N simulations
N = N.0*100
q.0 = 1/N

#number of simulations per inversion length
sims = 100*N
sim.results = rep(0, 10)

for(k in 1:10){
  x = k/10 - 0.05 #inversion length
  
  #number of generations per simulation (time until inversion is effectively neutral Ns_t = 1)
  gens = round(max(1, (1/s.d)*log(N*U*x)))
  
  #cumulative fixations
  sum.fix = 0
  
  for(j in 1:sims){
    q = q.0
    t = 0
    
    while((t < gens) & (q > 0)){
      #expected frequency after selection
      w.avg = q*exp(-U*x*(1 - exp(-s.d*t))) + (1 - q)*exp(-U*x)
      F.sel = q + q*(1 - q)*(exp(-U*x*(1 - exp(-s.d*t))) - exp(-U*x))/w.avg
      
      #binomial sampling
      q = rbinom(1, N, F.sel)/(N)
      t = t + 1
    }
    
    #sum fixations
    sum.fix = sum.fix + rbinom(1, 1, q)
  }
  
  sim.results[k] = exp(-U*x/s.d)*sum.fix/sims
  print(k)
}

large.N = sim.results

points((1:10/10 - 0.05),large.N*100, pch = 19)


#results using gamma DFE with shape k and k*theta = 0.02, and N = 10^4
#points((1:10/10 - 0.05), small.k, pch = 16, col = "orange")
#points((1:10/10 - 0.05),inter.k, pch = 16, col = "skyblue")
#points((1:10/10 - 0.05),large.k, pch = 16, col = "green")
