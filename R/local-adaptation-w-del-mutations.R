N = 10^6 #population size
U.d = 0.05 #deleterious mutation rate per chromosome arm
s.d = 0.01 #fitness effect of deleterious allele
m = 0.03 #migration rate
s.1 = 0.15 #selection at 1st locus
s.2 = 0.15 #selection at 2nd locus

sim.runs = 5*10^5 #simulations per inversion length

sim.fixations.nodel = rep(0, 10) #fixed inversions that were initially mutation-free
sim.fixations.wdel = rep(0, 10) #fixed inversions that carry mutations
sim.fixations.tot = rep(0, 10) #total fixed inversions
sim.lengths = rep(0, 10)

for(j in 1:10){
  x = j/10 - 0.05 #inversion length
  fixations.nodel = 0
  fixations.wdel = 0
  
  for(i in 1:sim.runs){
    r = x*abs(runif(1) - runif(1)) #assumes one crossover per chromosome arm per meiosis
    
    i.max = floor(-log(1 + m)/log(1 - s.d)) #maximum number of mutations an adaptive inversion can carry
    
    #deleterious mutations on the inversion
    del.mut = rpois(1, U.d*x/s.d)
    
    if(del.mut > i.max){
      fixations.nodel = fixations.nodel + 0
    }else{
      
      w.1 = 1 - s.1 - s.2
      w.2 = 1 - s.2
      w.3 = 1 - s.1
      w.4 = 1
      w.5 = (1 - s.d)^del.mut
      
      #Initial conditions of simulation (no LD approximation)
      x.1 = (m/s.1)*(m/s.2)
      x.2 = (1 - m/s.1)*m/s.2
      x.3 = (1 - m/s.2)*m/s.1
      x.4 = (1 - m/s.1)*(1 - m/s.2)
      x.5 = 0
      
      #deleterious allele frequencies and LD in population 1
      p.1 = x.1 + x.3
      p.2 = x.1 + x.2
      D = x.1*x.4 - x.2*x.3
      
      #deterministic migration-selection balance (100 generations burn-in)
      for(l in 1:100){
        #frequencies after migration
        x.1m = x.1*(1 - m) + m
        x.2m = x.2*(1 - m)
        x.3m = x.3*(1 - m)
        x.4m = x.4*(1 - m)
        x.5m = x.5*(1 - m)
        
        #deterministic frequencies after selection
        w.avg = x.1m*w.1 + x.2m*w.2 + x.3m*w.3 + x.4m*w.4 + x.5m*w.5
        x.1s = x.1m*w.1/w.avg
        x.2s = x.2m*w.2/w.avg
        x.3s = x.3m*w.3/w.avg
        x.4s = x.4m*w.4/w.avg
        x.5s = x.5m*w.5/w.avg
        
        #offspring of the next generation
        x.1 = x.1s - r*(x.1s*x.4s - x.2s*x.3s)
        x.2 = x.2s + r*(x.1s*x.4s - x.2s*x.3s)
        x.3 = x.3s + r*(x.1s*x.4s - x.2s*x.3s)
        x.4 = x.4s - r*(x.1s*x.4s - x.2s*x.3s)
        x.5 = x.5s
        
        #vectors of allele frequency and LD
        p.1 = c(p.1, x.1 + x.3)
        p.2 = c(p.2, x.1 + x.2)
        D = c(D, x.1*x.4 - x.2*x.3)
      }
      
      #initial conditions per simulation run
      x.1 = x.1*(1 - 1/N)
      x.2 = x.2*(1 - 1/N)
      x.3 = x.3*(1 - 1/N)
      x.4 = x.4*(1 - 1/N)
      x.5 = 1/N
      
      time = 0
      
      est.cutoff = 0.9
      
      while(x.5 > 0 & x.5 < est.cutoff){
        
        #frequencies after migration
        x.1m = x.1*(1 - m) + m
        x.2m = x.2*(1 - m)
        x.3m = x.3*(1 - m)
        x.4m = x.4*(1 - m)
        x.5m = x.5*(1 - m)
        
        #deterministic frequencies after selection
        w.avg = x.1m*w.1 + x.2m*w.2 + x.3m*w.3 + x.4m*w.4 + x.5m*w.5*exp(U.d*x*exp(-s.d*time))
        x.1s = x.1m*w.1/w.avg
        x.2s = x.2m*w.2/w.avg
        x.3s = x.3m*w.3/w.avg
        x.4s = x.4m*w.4/w.avg
        x.5s = x.5m*w.5*exp(U.d*x*exp(-s.d*time))/w.avg
        
        #genetic drift (sample N adults)
        x.drift = rmultinom(1, N, c(x.1s, x.2s, x.3s, x.4s, x.5s))/N
        x.1s = x.drift[1]
        x.2s = x.drift[2]
        x.3s = x.drift[3]
        x.4s = x.drift[4]
        x.5s = x.drift[5]
        
        #offspring of the next generation
        x.1 = x.1s - r*(x.1s*x.4s - x.2s*x.3s)
        x.2 = x.2s + r*(x.1s*x.4s - x.2s*x.3s)
        x.3 = x.3s + r*(x.1s*x.4s - x.2s*x.3s)
        x.4 = x.4s - r*(x.1s*x.4s - x.2s*x.3s)
        x.5 = x.5s
        
        time = time + 1
      }
      
      fixations.wdel = fixations.wdel + sign(x.5)*sign(del.mut)
      fixations.nodel = fixations.nodel + sign(x.5)*(1 - sign(del.mut))
    }
  }
  
  Pr.fix.nodel = fixations.nodel/sim.runs
  Pr.fix.wdel = fixations.wdel/sim.runs
  Pr.fix.tot = (fixations.nodel + fixations.wdel)/sim.runs
  print(x)
  
  sim.fixations.nodel[j] = Pr.fix.nodel
  sim.fixations.wdel[j] = Pr.fix.wdel
  sim.fixations.tot[j] = Pr.fix.tot
  sim.lengths[j] = x
}

par(mfrow = c(1, 1))
#probability that inversions capture the favoured genotype and become established (est. = reach a freq > 0.5)

#prediction inversions that are initially mutation-free
analytical.0del = function(x){
  x^2*2*(m*((2*(s.1 + s.2)*((s.1 + s.2) + x)*log((s.1 + s.2)/((s.1 + s.2) + x)))/x^2 + 2*(s.1 + s.2)/x + 1))*(1 + U.d*x/(1 - (1 - (m*((2*(s.1 + s.2)*((s.1 + s.2) + x)*log((s.1 + s.2)/((s.1 + s.2) + x)))/x^2 + 2*(s.1 + s.2)/x + 1)))*exp(-s.d)))*exp(-U.d*x/s.d)
}

#prediction with U.d = 0
analytical.0del.limit = function(x){
  x^2*2*(m*((2*(s.1 + s.2)*((s.1 + s.2) + x)*log((s.1 + s.2)/((s.1 + s.2) + x)))/x^2 + 2*(s.1 + s.2)/x + 1))
}

curve(analytical.0del.limit, 0, 1, ylim = c(10^-5, 10^-1), xlab = "inversion length (x)", ylab = "fixation probability", log = "y")
points(sim.lengths, sim.lengths^2*sim.fixations.nodel, pch = 16, col = "grey")
curve(analytical.0del, add = TRUE, col = "GREY", lwd = 3)
points(sim.lengths, sim.lengths^2*sim.fixations.tot, pch = 16)
