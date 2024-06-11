using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(9999)

function main2(N,jl,operator::String)

	sites = siteinds("S=1/2", N)
	
	os = getOperator(operator, N)
	
	operator = MPO(os, sites)
	
	psi0 = randomMPS(sites; linkdims=10)
	
	nsweeps = 10
	maxdim = [10,20,100,100,200]
	
	cutoff = [1e-6]
	
	energy, psi = dmrg(operator, psi0; nsweeps, maxdim, cutoff)
	
	return [energy, psi]
end

#TODO better way of defining this, rework entirely
function getOperator(op::String, nsites::Integer)
	op = lowercase(op)		#force case insensitive
	
	if (cmp(op,"hamiltonian")== 0 || cmp(op,"h")== 0 || cmp(op,"ham")== 0)
	#this is a hamiltonian
		os = OpSum()
		for j in 1:(nsites-1)
			os += jl,"Sz", j, "Sz", j+1
			os += jl*0.5, "S+", j, "S-", j+1
			os += jl*0.5, "S-", j, "S+", j+1
		end
		
	elseif (occursin("spin",op))
	#this is a spin expectation on site N
	
	#now get the int here IG taking format spin5, spin(5), or spins(a,b,c,d)
	#using ~regex~
	#there must be a better way
	
	#SzM = op("Sz", sites, M)
	
	end
	
	return ArgumentError("operator specified was not in the list of known operators")
end

main2(2,1,"H")
