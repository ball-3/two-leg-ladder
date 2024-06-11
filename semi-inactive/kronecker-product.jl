function kron(matrices...)

	len = length(matrices)

	if len == 0
		return nothing
	elseif len == 1
		return matrices[1]
	end
	
	rowsA = size(matrices[1],1)
	columnsA = size(matrices[1],2)
	matrixA = matrices[1]
	
	for i = 2:(len)
		rowsB = size(matrices[i],1)
		columnsB = size(matrices[i],2)
		matrixB = matrices[i]
		
		rowsAB = rowsA*rowsB
		colsAB = columnsA*columnsB
		matrixAB = zeros(Complex{Float64},rowsAB,colsAB)
		
		cAB = 0
		for cA = 1:columnsA
			for cB = 1:columnsB
				cAB += 1
				rAB = 0
				for rA = 1:rowsA
					for rB = 1:rowsB
						rAB += 1
						matrixAB[rAB,cAB] = matrixA[rA,cA]*matrixB[rB,cB]
					end
				end
			end
		end
		
		rowsA = rowsAB
		columnsA = colsAB
		matrixA = matrixAB
	end
	return matrixA
end
