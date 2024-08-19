
using LinearAlgebra
using Printf

function ED(H)
	eigs = eigvals(H)
	n = length(eigs)
	vecM = eigvecs(H)
	vecs = [zeros(n) for i in 1:n]
	
	for i in 1:n
		vecs[i] = vecM[:,i]
	end

	return [eigs,vecs]
end

function displayED(H, precision)
	eigs, vecs = ED(H)
	n = length(eigs)
	
	for i in 1:n
		@printf("Eigenpair %d of %d:\n  %.*f [%.*f", i, n, precision, eigs[i], precision, vecs[i][1])
		len = length(eigs[i])
		
		for j in 2:n
			#@printf("\n   %*.5f", (5+len), vecs[i][j])
			@printf("\n   %*.*f", 2*precision+5, precision,vecs[i][j])
		end
		@printf("]\n\n")
	end
end

function gsE(H)
	return minimum(eigvals(H))
end

function gsE(H,psi)
	vals = ED(H)
	eigs = Vector{ComplexF64}(vals[1])
	vecs = vals[2]
	n = length(psi)
	c = [psi'*vecs[i] for i in 1:n]
	E = 0
	#psiApprox = [0 + 0*im for i in 1:n]
	for i in 1:n
		E += c[i]*conj(c[i])*eigs[i]
		#psiApprox += c[i]*vecs[i]
	end
	#display(psi)
	#println("is exact")
	#display(psiApprox)
	#println("is approximate")

	return E
end