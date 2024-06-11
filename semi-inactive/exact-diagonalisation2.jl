include("mpo.jl")
include("kronecker-product.jl")

#this function generates the visual representation, not the matrix math representation
#=function ket(i,n)	#up to 64 particles but youre not going to do that anyways
	N = 2^n
	s = bitstring(Int64(i-1))
	s = last(s,n)
	S = split(s,"")
	ket = zeros(n)
	for j = 1:n
		if S[j] == "1"
			ket[j] = 1
		end
	end
	return ket
end=#

function ket(i,N)
	ket = zeros(N)
	ket[i] = 1
	return ket	
end

function bra(ket::Vector)
	return ket'
end

function ED(n::Integer)
	N = 2^n
	
	H = MPO(n)
	display(H)
	kets = [ket(i,N) for i in 1:N]
	bras = [bra(kets[i]) for i in 1:N]
	
	A = zeros(N,N)
	
	for i = 1:N
		for j = 1:N
			A[i,j] = bras[i]*H*kets[j]
		end
	end
	return A
end
