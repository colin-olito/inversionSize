#error function
erf = function(y){2*pnorm(y*sqrt(2))-1}

#parameters
U = 0.02
N = 10^4
s.d = 0.01
s = 0.0001
2*N*s

#number of simulation runs per parameter set and x
sim.runs = 1000*N

#initial inversion frequency
q.0 = 1/(2*N)

#predicted fixation probability
prob.fix = function(size){
  (erf(sqrt(N*s*size)) - erf(sqrt(N*s*size)*(1 - 2*q.0)))/(2*erf(sqrt(N*s*size)))
}

#simulations with small load
sizes = 1:10/10 - 0.05
simulated = rep(0, length(sizes))

for(j in 1:length(sizes)){
  #inversion size
  x = sizes[j]
  
  #fixation events (initially 0)
  successes = 0
  
  for(i in 1:sim.runs){
    t = 0
    q = q.0
    
    while(q*(1 - q) > 0){
      #expected frequency after selection
      w.avg = q^2*exp(-2*U*x*(1 - exp(-s.d*t))) + 2*q*(1 - q)*(1 - s*x)*exp(-U*x*(2 - exp(-s.d*t))) + (1 - q)^2*exp(-2*U*x)
      F.sel = q + q*(1 - q)*(q*exp(-2*U*x*(1 - exp(-s.d*t))) + (1 - 2*q)*(1 - s*x)*exp(-U*x*(2 - exp(-s.d*t))) - (1 - q)*exp(-2*U*x))/w.avg
      
      #binomial sampling
      q = rbinom(1, 2*N, F.sel)/(2*N)
      t = t + 1
    }
    
    if(q > 0){successes = successes + 1} else{successes = successes + 0}
  }
  
  print(c(x, s, U, prob.fix(x), successes/sim.runs))
  simulated[j] = successes/sim.runs
}
Pr.mut.free = exp(-U*sizes/s.d)
sim.w.load.1 = Pr.mut.free*simulated

y = 1:1000/1000
par(mfrow = c(1, 1))
plot(y, prob.fix(y), type = "l", xlab = "inversion length (x)", ylab = "fixation probability", ylim = c(0, 1.25/(2*N)), xlim = c(0, 1), lwd = 3)
points(sizes,sim.w.load.1, pch = 21, bg = "WHITE")

#simulations with intermediate load
U = 0.05
simulated = rep(0, length(sizes))
for(j in 1:length(sizes)){
  #inversion size
  x = sizes[j]
  
  #fixation events (initially 0)
  successes = 0
  
  for(i in 1:sim.runs){
    t = 0
    q = q.0
    
    while(q*(1 - q) > 0){
      #expected frequency after selection
      w.avg = q^2*exp(-2*U*x*(1 - exp(-s.d*t))) + 2*q*(1 - q)*(1 - s*x)*exp(-U*x*(2 - exp(-s.d*t))) + (1 - q)^2*exp(-2*U*x)
      F.sel = q + q*(1 - q)*(q*exp(-2*U*x*(1 - exp(-s.d*t))) + (1 - 2*q)*(1 - s*x)*exp(-U*x*(2 - exp(-s.d*t))) - (1 - q)*exp(-2*U*x))/w.avg
      
      #binomial sampling
      q = rbinom(1, 2*N, F.sel)/(2*N)
      t = t + 1
    }
    
    if(q > 0){successes = successes + 1} else{successes = successes + 0}
  }
  
  print(c(x, s, U, prob.fix(x), successes/sim.runs))
  simulated[j] = successes/sim.runs
}
Pr.mut.free = exp(-U*sizes/s.d)
sim.w.load.2 = Pr.mut.free*simulated

points(sizes,sim.w.load.2, pch = 19, col = "GREY")


#simulations with large load
U = 0.1
simulated = rep(0, length(sizes))
for(j in 1:length(sizes)){
  #inversion size
  x = sizes[j]
  
  #fixation events (initially 0)
  successes = 0
  
  for(i in 1:sim.runs){
    t = 0
    q = q.0
    
    while(q*(1 - q) > 0){
      #expected frequency after selection
      w.avg = q^2*exp(-2*U*x*(1 - exp(-s.d*t))) + 2*q*(1 - q)*(1 - s*x)*exp(-U*x*(2 - exp(-s.d*t))) + (1 - q)^2*exp(-2*U*x)
      F.sel = q + q*(1 - q)*(q*exp(-2*U*x*(1 - exp(-s.d*t))) + (1 - 2*q)*(1 - s*x)*exp(-U*x*(2 - exp(-s.d*t))) - (1 - q)*exp(-2*U*x))/w.avg
      
      #binomial sampling
      q = rbinom(1, 2*N, F.sel)/(2*N)
      t = t + 1
    }
    
    if(q > 0){successes = successes + 1} else{successes = successes + 0}
  }
  
  print(c(x, s, U, prob.fix(x), successes/sim.runs))
  simulated[j] = successes/sim.runs
}

Pr.mut.free = exp(-U*sizes/s.d)
sim.w.load.3 = Pr.mut.free*simulated

points(sizes,sim.w.load.3, pch = 19)

