module Hamiltonians

include("kronecker-product.jl")
include("pauli-matrices.jl")

export heisenbergOneHalf, ladderOneHalf

function heisenbergOneHalf(spins::Integer, jl)
#a hamiltonian of H = JΣᵢSi(dot)Si+1 = JΣᵢ(1/2*SxiSxi+1 + 1/2*SyiSyi+1 SziSz+1) and all particles spin 1/2
	dim = 2^spins
	H = zeros(dim,dim)
	hbar = 1
	j = (spins-2)

	for i = 0:(spins-2)
		H += jl*kron(nthI(i),(hbar/2)*σx,(hbar/2)*σx,nthI(j))
		H += jl*kron(nthI(i),(hbar/2)*σy,(hbar/2)*σy,nthI(j))
		H += jl*kron(nthI(i),(hbar/2)*σz,(hbar/2)*σz,nthI(j))
		j -= 1
	end
	return H
end

function ladderOneHalf(pairs::Integer, Jl, Jr)
	#settings:
	hbar = 1
	
	N = 2*pairs
	dim = 2^(N)
	H = zeros(Complex{Float64},dim,dim)

	j = N-2
	for i = 0:2:(N-2)
		#intra dimer interaction (Jr)
		H += Jr*kron(nthI(i),(hbar/2)*σx,(hbar/2)*σx,nthI(j))
		H += Jr*kron(nthI(i),(hbar/2)*σy,(hbar/2)*σy,nthI(j))
		H += Jr*kron(nthI(i),(hbar/2)*σz,(hbar/2)*σz,nthI(j))
		j -= 2
	end

	j = N-4
	for i = 0:2:(N-4)
		#inter dimer interaction (Jl)
			H += Jl*kron(nthI(i),nthI(1),(hbar/2)*σx,nthI(1),(hbar/2)*σx,nthI(j))
			H += Jl*kron(nthI(i),nthI(1),(hbar/2)*σy,nthI(1),(hbar/2)*σy,nthI(j))
			H += Jl*kron(nthI(i),nthI(1),(hbar/2)*σz,nthI(1),(hbar/2)*σz,nthI(j))
			
			H += Jl*kron(nthI(i),(hbar/2)*σx,nthI(1),(hbar/2)*σx,nthI(1),nthI(j))
			H += Jl*kron(nthI(i),(hbar/2)*σy,nthI(1),(hbar/2)*σy,nthI(1),nthI(j))
			H += Jl*kron(nthI(i),(hbar/2)*σz,nthI(1),(hbar/2)*σz,nthI(1),nthI(j))
		j -= 2
	end
	return H
end

function nthI(times::Integer)		#surely there is a better way to name this function
	if times == 0 return 1 end
	working = I2
	
	for i = 1:(times-1)
		working = kron(working,I2)
	end
	return convert(Array{ComplexF64},working)
end

end
