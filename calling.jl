include("./utility-functions/comparison-functions.jl")
include("plotting.jl")

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

function callingABCTest(nPairs, ED::Bool, MPS::Bool, Jl, Jr)
	sigma = 0.5*[σx,σy,σz]
	sigmaStr = ["Sx","Sy","Sz"]
	i = 1

	expA, expB, expC = 0, 0, 0

	if(ED)
		groundState = ED(Hamiltonians.ladderOneHalf(nPairs,Jl,Jr))[2][1]

		for j in 1:3
			op = [sigma[j],sigma[j]]
			opStr = "<"*sigmaStr[j]*sigmaStr[j]*">"

			expA = CustOp(groundState,op,[i,i+1])
			expB = CustOp(groundState,op,[i,i+2])
			expC = CustOp(groundState,op,[i+1,i+3])

			exp = expA - 0.5*expB - 0.5*expC
			#@printf("Using ED for operator %s  : %f = %f - 1/2 %f - 1/2 %f\n", opStr, exp, expA, expB, expC)

		end
	end

	if(MPS)
		groundState = dmrgLadder(nPairs, Jl, Jr)
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

			#TODO this is a temp change to fit my plot fct, edit later?
			return [expA,expB,expC]
			#exp = expA - 0.5*expB - 0.5*expC
			#@printf("Using MPS for operator %s : %f = %f - 1/2 %f - 1/2 %f\n", opStr, exp, expA, expB, expC)
		end
	end
end

function makeABCPlot(numDataSets,maxPairs)
	title = "Spin-Spin correlations as a function of Coupling Strength"
	xTitle = "number of site pairs"
	yTitle = "spin-spin correlation"
	save = true
	
	labels = ["$i" for i in 1:2*numDataSets]
	multiDataAy = [[0.0 for i in 1:maxPairs] for j in 1:numDataSets]
	multiDataBy = [[0.0 for i in 1:maxPairs] for j in 1:numDataSets]
	multiDataCy = [[0.0 for i in 1:maxPairs] for j in 1:numDataSets]
	datax = [i for i in 1:maxPairs]

	#TODO initialise data
	Jr = 1
	Jl = 1000
	for k in 1:numDataSets

		Jl = Jl/(10)
		JrbyJl = Jr/Jl
		labels[2*k-1] = "Jr/Jl : $(JrbyJl)"

		dataAy = [0.0 for i in 1:maxPairs]
		dataBy = [0.0 for i in 1:maxPairs]
		dataCy = [0.0 for i in 1:maxPairs]

		for j in 1:maxPairs
			display(j)
			display(dataAy)
			dataAy[j], dataBy[j], dataCy[j] = callingABCTest(j,false,true,Jl,Jr)
		end

		multiDataAy[k] = dataAy
		multiDataBy[k] = dataBy
		multiDataCy[k] = dataCy
	end
	myScatterPlot(title*"A",true,xTitle,yTitle,labels,datax,multiDataAy[1],datax,multiDataAy[2],datax,multiDataAy[3],datax,multiDataAy[4])
	myScatterPlot(title*"B",true,xTitle,yTitle,labels,datax,multiDataBy[1],datax,multiDataBy[2],datax,multiDataBy[3],datax,multiDataBy[4])
	myScatterPlot(title*"C",true,xTitle,yTitle,labels,datax,multiDataCy[1],datax,multiDataCy[2],datax,multiDataCy[3],datax,multiDataCy[4])
end


