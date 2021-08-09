N = 10^6 #population size
U.d = 0.1 #deleterious mutation rate per chromosome arm
s.het = 0.002 #heterozygous fitness effect of deleterious allele
h = 0.25 #dominance coefficient
s.d = s.het/h #homozygous fitness effect of deleterious allele
s.b = 20*s.het #homozygous fitness benefit of the inversion
U.d/s.het

i.max = floor(-log(1 + s.b/2)/log(1 - s.d*h)) #maximum number of mutations an adaptive inversion can carry
i.bal = floor(log(1 - s.b/(2 + 2*s.b))/log(1 - s.d*(1 - h)/(1 - s.d*h))) #threshold for balancing selection
c(i.max, i.bal)

#simulations per inversion length
if(s.b > 0){sim.runs = 10^5}else{sim.runs = 200*N}

sim.fixations = rep(0, 10) #fixed inversions
sim.balanced = rep(0, 10) #balanced inversions
sim.lengths = rep(0, 10)

for(j in 1:10){
  x = j/10 - 0.05 #inversion length
  fixations = 0
  balanced = 0
  
  for(i in 1:sim.runs){
    #deleterious mutations on the inversion
    del.mut = rpois(1, U.d*x/s.het)
    
    #fitness per genotype
    w.SS = 1
    w.IS = (1 + s.b/2)*(1 - s.het)^del.mut
    w.II = (1 + s.b)*(1 - s.d)^del.mut
    
    #starting frequencies at birth
    p = 1 - 1/(2*N)
    q = 1/(2*N)
    
    time = 0
    
    if(del.mut <= i.max){
      q.eq = (w.IS - w.SS)/(2*w.IS - w.II - w.SS) #intermediate equilibrium (if there is one)
      
      if(0 < q.eq & q.eq < 1){
        
        while(q > 0 & q < q.eq){
          w.avg = p^2*w.SS + 2*p*q*w.IS*exp(U.d*x*exp(-s.het*time)) + q^2*w.II*exp(2*U.d*x*exp(-s.het*time))
          SS.sel = p^2*w.SS/w.avg
          IS.sel = 2*p*q*w.IS*exp(U.d*x*exp(-s.het*time))/w.avg
          II.sel = q^2*w.II*exp(2*U.d*x*exp(-s.het*time))/w.avg
          
          drift = rmultinom(1, N, c(SS.sel, IS.sel, II.sel))/N
          
          p = drift[1] + drift[2]/2
          q = drift[2]/2 + drift[3]
          
          time = time + 1
        }
        
        balanced = balanced + sign(q)
      } else{
        
        while(q*(1 - q) > 0){
          w.avg = p^2*w.SS + 2*p*q*w.IS*exp(U.d*x*exp(-s.het*time)) + q^2*w.II*exp(2*U.d*x*exp(-s.het*time))
          SS.sel = p^2*w.SS/w.avg
          IS.sel = 2*p*q*w.IS*exp(U.d*x*exp(-s.het*time))/w.avg
          II.sel = q^2*w.II*exp(2*U.d*x*exp(-s.het*time))/w.avg
          
          drift = rmultinom(1, N, c(SS.sel, IS.sel, II.sel))/N
          
          p = drift[1] + drift[2]/2
          q = drift[2]/2 + drift[3]
          
          time = time + 1
        }
        
        fixations = fixations + q
      }
      
    } else{del.mut = del.mut}
    
  }
  
  Pr.fixation = fixations/sim.runs
  Pr.balancing = balanced/sim.runs
  print(x)
  
  sim.fixations[j] = Pr.fixation
  sim.balanced[j] = Pr.balancing
  sim.lengths[j] = x
  
}

par(mfrow = c(1, 1))

if(s.b > 0){
  #Peischel & Kirkpatrick meets Manning and Thompson
  lengths = rep(0, 1001)

  Predicted.pr.est = rep(0, 1001)
  Predicted.pr.fix = rep(0, 1001)
  Predicted.pr.bal = rep(0, 1001)
  
  for(j in 1:1001){
    lengths[j] = (j - 1)/1000
    x = (j - 1)/1000
    i.mut = 0
    PKMT.est = 0
    
    while(i.mut <= i.max){
      s.i = (1 + s.b/2)*(1 - s.het)^i.mut - 1
      PKMT.est = PKMT.est + 2*(s.i)*(1 + U.d*x/(1 - (1 - s.i)*exp(-s.het)))*exp(-U.d*x/s.het)*(U.d*x/s.het)^i.mut/factorial(i.mut)
      i.mut = i.mut + 1
    }
    Predicted.pr.est[j] = PKMT.est
  }
  
  for(j in 1:1001){
    lengths[j] = (j - 1)/1000
    x = (j - 1)/1000
    i.mut = 0
    PKMT.fix = 0
    while(i.mut <= i.bal){
      s.i = (1 + s.b/2)*(1 - s.het)^i.mut - 1
      PKMT.fix = PKMT.fix + 2*(s.i)*(1 + U.d*x/(1 - (1 - s.i)*exp(-s.het)))*exp(-U.d*x/s.het)*(U.d*x/s.het)^i.mut/factorial(i.mut)
      i.mut = i.mut + 1
    }
    Predicted.pr.fix[j] = PKMT.fix
  }  
  
  for(j in 1:1001){
    lengths[j] = (j - 1)/1000
    x = (j - 1)/1000
    i.mut = i.bal + 1
    PKMT.bal = 0
    while(i.mut <= i.max){
      s.i = (1 + s.b/2)*(1 - s.het)^i.mut - 1
      PKMT.bal = PKMT.bal + 2*(s.i)*(1 + U.d*x/(1 - (1 - s.i)*exp(-s.het)))*exp(-U.d*x/s.het)*(U.d*x/s.het)^i.mut/factorial(i.mut)
      i.mut = i.mut + 1
    }
    Predicted.pr.bal[j] = PKMT.bal  
  }
  
}else{
  lengths = 0:1000/1000
  Predicted.pr.est = rep(1/(2*N), 1001)
}

plot(lengths, Predicted.pr.est, type = "l", xlab = "inversion length (x)", ylab = "probability", ylim = c(0, max(1.2/(2*N), s.b)), lty = 3)
lines(lengths, Predicted.pr.bal, col = "RED")
lines(lengths, Predicted.pr.fix)
points(sim.lengths, sim.balanced, pch = 19, col = "RED")
points(sim.lengths, sim.fixations, pch = 19)
points(sim.lengths, sim.balanced + sim.fixations, lwd = 2)
