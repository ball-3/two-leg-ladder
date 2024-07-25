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
	sum = 0
	tn = 2*n
	for i in 1:(tn-1)
		#sum += CustOp(groundState,[S⁺*S⁻, S⁻*S⁺],[1,2])
	end
	display(sum/tn)


	groundState = dmrgLadder(n, 1, 1)
	groundStateWFCT = groundState[2]
	sites = groundState[3]

	display(CustOp(sites,groundStateWFCT,["S- * S+"],[1]))
end

callingABCTest(1)