N = 10^6 #population size
U.d = 0.1 #deleterious mutation rate per chromosome arm
s.d = 0.01 #fitness effect of deleterious allele
s.b = s.d/2 #benefit of the inversion

i.max = floor(-log(1 + s.b)/log(1 - s.d)) #maximum number of mutations an adaptive inversion can carry

#simulations per inversion length
if(s.b > 0){sim.runs = 10^5}else{sim.runs = 100*N}

sim.fixations.nodel = rep(0, 10) #fixed inversions that were initially mutation-free
sim.fixations.wdel = rep(0, 10) #fixed inversions that carry mutations
sim.fixations.tot = rep(0, 10) #total fixed inversions
sim.lengths = rep(0, 10)

for(j in 1:10){
  x = j/10 - 0.05 #inversion length
  fixations.nodel = 0
  fixations.wdel = 0
  
  for(i in 1:sim.runs){
    #deleterious mutations on the inversion
    del.mut = rpois(1, U.d*x/s.d)
    
    #fitness per genotype
    w.1 = 1
    w.2 = (1 + s.b)*(1 - s.d)^del.mut
    
    #starting frequencies at birth
    p = 1 - 1/N
    q = 1/N
    
    time = 0
    
    while(q*(1 - q) > 0){
      w.avg = p*w.1 + q*w.2*exp(U.d*x*exp(-s.d*time))
      p.sel = p*w.1/w.avg
      q.sel = q*w.2*exp(U.d*x*exp(-s.d*time))/w.avg
      
      drift = rmultinom(1, N, c(p.sel, q.sel))/N
      
      p = drift[1]
      q = drift[2]
      
      time = time + 1
    }
    
    fixations.wdel = fixations.wdel + q*sign(del.mut)
    fixations.nodel = fixations.nodel + q*(1 - sign(del.mut))
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

if(s.b > 0){
  #probability of inversion becoming fixed
  analytical.0del = function(x){
    2*s.b*(1 + U.d*x/(1 - (1 - s.b)*exp(-s.d)))*exp(-U.d*x/s.d)
  }
  
  #Peischel & Kirkpatrick meets Manning and Thompson
  
  lengths = rep(0, 51)
  
  Predicted.pfix = rep(0, 51)
  for(j in 1:51){
    lengths[j] = (j - 1)/50
    x = (j - 1)/50
    i.mut = 0
    PKMT = 0
    while(i.mut <= i.max){
      s.i = (1 + s.b)*(1 - s.d)^i.mut - 1
      PKMT = PKMT + 2*(s.i)*(1 + U.d*x/(1 - (1 - s.i)*exp(-s.d)))*exp(-U.d*x/s.d)*(U.d*x/s.d)^i.mut/factorial(i.mut)
      i.mut = i.mut + 1
    }
    Predicted.pfix[j] = PKMT
  }
}else{
  lengths = 0:50/50
  Predicted.pfix = rep(1/N, 51)
  
  #probability of inversion becoming fixed
  analytical.0del = function(x){
    1/N + 0*x
  }
  
}

plot(lengths, Predicted.pfix, type = "l", xlab = "inversion length (x)", ylab = "fixation probability", ylim = c(0, max(1.2/N, 2*s.b)))
curve(analytical.0del, add = TRUE, col = "grey")
points(sim.lengths, sim.fixations.nodel, pch = 16, col = "grey")
points(sim.lengths, sim.fixations.tot, pch = 16)
