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
