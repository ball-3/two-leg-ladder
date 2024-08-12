using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(6969)

function dmrgLadder(numPairs, Jl, Jr)

	N = numPairs*2#TODO added this july 2, why whasnt it here?
	
	sites = siteinds("S=1/2", N)
	
	os = OpSum()
	for j in 1:(N-2)
		#inter dimer
		os += Jl, "Sz", j, "Sz", j+2
		os += Jl*0.5, "S+", j, "S-", j+2
		os += Jl*0.5, "S-", j, "S+", j+2
	end
	
	for j in 1:2:(N-1)	
		#intra dimer
		os += Jr, "Sz", j, "Sz", j+1
		os += Jr*0.5, "S+", j, "S-", j+1
		os += Jr*0.5, "S-", j, "S+", j+1
	end
	
	op = MPO(os, sites)
	
	psi0 = randomMPS(sites; linkdims=200)#changed from 10
	
	nsweeps = 20#changed from 10
	maxdim = [200]#increased
	
	cutoff = [1e-12]#changed from 1e-6
	
	energy, psi = dmrg(op, psi0; nsweeps, maxdim, cutoff, outputlevel = 0)
	
	#@printf("Final energy = %.12f, n = %d, Jl = %f, Jr = %f\n", energy, N, Jl, Jr)
	
	return [energy, psi, sites]
end

function timeEvolutionChain(numSites, Jl, Jr)
  N = numSites
  cutoff = 1E-8
  tau = 0.1
  ttotal = 5.0

  # Make an array of 'site' indices
  s = siteinds("S=1/2", N; conserve_qns=true)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ITensor[]
  for j in 1:2:(N - 3)
    s1 = s[j]
    s2 = s[j + 1]
    s3 = s[j + 2]
    s4 = s[j + 3]
    hj =
      Jr / 2 * (op("S-", s1) * op("S+", s2) + op("S+", s1) * op("S-", s2) + 
                op("S-", s3) * op("S+", s4) + op("S+", s3) * op("S-", s4) +
                op("Sz", s1) * op("Sz", s2) + op("Sz", s3) * op("Sz", s4))+
      Jl / 2 * (op("S-", s1) * op("S+", s3) + op("S+", s1) * op("S-", s3) + 
                op("S-", s2) * op("S+", s4) + op("S+", s2) * op("S-", s4) +
                op("Sz", s1) * op("Sz", s3) + op("Sz", s2) * op("Sz", s4))
      #op("Sz", s1) * op("Sz", s2) +
      #1 / 2 * op("S+", s1) * op("S-", s2) +
      #1 / 2 * op("S-", s1) * op("S+", s2)
    Gj = exp(-im * tau / 2 * hj)
    push!(gates, Gj)
  end
  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates, reverse(gates))

  # Initialize psi to be a product state (alternating up and down)
  psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")

  c = div(N, 2) # center site

  # Compute and print <Sz> at each time step
  # then apply the gates to go to the next time
  for t in 0.0:tau:ttotal
    Sz = expect(psi, "Sz"; sites=c)
    println("$t $Sz")

    t≈ttotal && break

    psi = apply(gates, psi; cutoff)
    normalize!(psi)
  end

  return
end

timeEvolutionChain(10,1,1)