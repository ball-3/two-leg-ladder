include("./utility-functions/comparison-functions.jl")

include("./dmrg/heisenberg.jl")
include("./dmrg/ladder.jl")
include("./dmrg/operators.jl")

include("./exact-diagonalisation/Hamiltonians.jl")
include("./exact-diagonalisation/expectations.jl")
include("./exact-diagonalisation/exact-diagonalisation.jl")

#displayED(Hamiltonians.ladderOneHalf(2,1,1),5)

function callingNormalisationTest()
	n = 4
	m = 2^(2*n)
	
	for i in 1:m
		display(normalisationTest(ED(Hamiltonians.ladderOneHalf(n,1,1))[2][i]))
	end
end


function callingSpSmTest(n)
	groundState = ED(Hamiltonians.ladderOneHalf(n,1,1))[2][1]
	print("\nED expectation ")
	sum = 0
	tn = 2*n
	
	#for i in 1:(tn-1)
		#sum += CustOp(groundState,[S⁺*S⁻, S⁻*S⁺],[1,2])
		sum += CustOp(groundState,[S⁻*S⁺],[1])
	#end
	display(sum)

	
	groundState = dmrgLadder(n, 1, 1)
	groundStateWFCT = groundState[2]
	sites = groundState[3]
	print("\nDMRG expectation ")
	display(CustOp(sites,groundStateWFCT,["S- * S+"],[1]))
	#display(CustOp(sites,groundStateWFCT,["S+ * S-","S- * S+"],[1,2]))
	
end

function callingABCTest(n)
	groundState = ED(Hamiltonians.ladderOneHalf(n,1,1))[2][1]
	tn = 2*n
	sigma = [σx,σy,σz]


	for j in 1:3
		op = [sigma[j],sigma[j]]
		expA, expB, expC = 0, 0, 0

		for i in 1:2:(tn-1)
			expA += 2*CustOp(groundState,op,[i,i+1])/(tn-1)
		end
		for i in 1:(tn-3)
			expB += CustOp(groundState,op,[i,i+2])/(tn-3)
			expC += CustOp(groundState,op,[i,i+3])/(tn-3)
		end

		exp = expA - 0.5*expB - 0.5*expC
		print("ED for spin-spin $(j)\n")
		display(exp)

	end
	


	sigmaStr = ["Sx","Sy","Sz"]

	groundState = dmrgLadder(n, 1, 1)
	groundStateWFCT = groundState[2]
	sites = groundState[3]

	for j in 1:3
		op = [sigmaStr[j],sigmaStr[j]]
		expA, expB, expC = 0, 0, 0

		for i in 1:2:(tn-1)
			expA += 2*CustOp(sites,groundStateWFCT,op,[i,i+1],2)/(tn-1)
		end
		for i in 1:(tn-3)
			expB += CustOp(sites,groundStateWFCT,op,[i,i+2],2)/(tn-3)
			expC += CustOp(sites,groundStateWFCT,op,[i,i+3],2)/(tn-3)
		end

		exp = expA - 0.5*expB - 0.5*expC
		print("DMRG for spin-spin $(j)\n")
		display(exp)
	end
end

callingABCTest(1)
callingABCTest(2)