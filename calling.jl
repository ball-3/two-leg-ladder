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
		sum += CustOp(groundState,[0.5*σz,0.5*σz],[1,2])
	#end
	display(sum)

	
	groundState = dmrgLadder(n, 1, 1)
	groundStateWFCT = groundState[2]
	sites = groundState[3]
	print("\nDMRG expectation ")
	display(SzSz(sites,groundStateWFCT))
	#display(CustOp(sites,groundStateWFCT,["S+ * S-","S- * S+"],[1,2]))
	
end

function callingABCTest(nPairs)
	sigma = 0.5*[σx,σy,σz]
	sigmaStr = ["Sx","Sy","Sz"]
	i = 1

	expA, expB, expC = 0, 0, 0

	groundState = ED(Hamiltonians.ladderOneHalf(nPairs,1,1))[2][1]

	for j in 1:3
		op = [sigma[j],sigma[j]]
		opStr = "<"*sigmaStr[j]*sigmaStr[j]*">"

		expA = CustOp(groundState,op,[i,i+1])
		expB = CustOp(groundState,op,[i,i+2])
		expC = CustOp(groundState,op,[i+1,i+3])

		exp = expA - 0.5*expB - 0.5*expC
		@printf("Using ED for operator %s  : %f = %f - 1/2 %f - 1/2 %f\n", opStr, exp, expA, expB, expC)

	end

	groundState = dmrgLadder(nPairs, 1, 1)
	groundStateWFCT = groundState[2]
	sites = groundState[3]

	for j in 1:3
		op = [sigmaStr[j],sigmaStr[j]]
		opStr = "<"*sigmaStr[j]*sigmaStr[j]*">"

		expA = CustOp(sites,groundStateWFCT,op,[i,i+1],2)

		if (nPairs > 1)
			expB = CustOp(sites,groundStateWFCT,op,[i,i+2],2)
			expC = CustOp(sites,groundStateWFCT,op,[i+1,i+3],2)
		end

		exp = expA - 0.5*expB - 0.5*expC
		@printf("Using MPS for operator %s : %f = %f - 1/2 %f - 1/2 %f\n", opStr, exp, expA, expB, expC)
	end
end

callingABCTest(1)
callingABCTest(2)
callingABCTest(3)