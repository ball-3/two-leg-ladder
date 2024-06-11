using LinearAlgebra

function ED(H)	#rename
	eigs = eigvals(H)
	display(eigs)
	n = length(eigs)
	vecs = [zeros(n) for i in 1:n]
	
	for i = 1:n
		vecs[i][i] = eigs[i]
	end
	return [eigs,vecs]
end

function gsE(H)
	return minimum(eigvals(H))
end
