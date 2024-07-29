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

#function to plot figure 4 from ground state magnetic properties of spin ladder shaped quantum nanomagnet: exact diagonalisation study
#using my mpo methods to test accuracy
function MMMfig4cmp()
	J1 = 1
	J2 = 3
	yVals1 = MMMfig4helper(J1,J2)
	J1 = -1
	yVals2 = MMMfig4helper(J1,J2)

	title = "MMM-fig4-using-MPO"
	savePlot = true
	xTitle = "J2 / |J1|"
	yTitle = "E0 / |J1|"
	labels = ["J1 < 0";"J1 < 0";"J1 > 0";"J1 > 0"]
	xVals = [-3;-2;-1;0;1;2;3]

	myScatterPlot(title, savePlot, xTitle, yTitle, labels, xVals, yVals1, xVals, yVals2,[-3;-3],[-16;0])
end

function MMMfig4helper(Jl, Jr)
	result = [0.0 for i in 1:6]# result = ground state energy by absolute value of Jl, y values in fig 4
	for i in 1:6
		absJ1 = abs(Jl)
		result[i] = dmrgLadder(6, Jl, Jr)[1]/absJ1	# = E0/|J1|
		Jr -= 1
	end
	return result
end

function MMMfig5cmp()
	J1 = 1
	J2 = 4
	c12 = MMMfig5helper(J1,J2,1,2)
	c34 = MMMfig5helper(J1,J2,3,4)
	c56 = MMMfig5helper(J1,J2,5,6)
	c13 = MMMfig5helper(J1,J2,1,3)
	c35 = MMMfig5helper(J1,J2,3,5)
	c57 = MMMfig5helper(J1,J2,5,7)
	c14 = MMMfig5helper(J1,J2,1,4)
	c36 = MMMfig5helper(J1,J2,3,6)
	c58 = MMMfig5helper(J1,J2,5,8)

	title = "MMM-fig5-using-MPO-A"
	savePlot = true
	xTitle = "J2 / |J1|"
	yTitle = "correlation"
	labels = ["c12";"";"c34";"";"c56";"";"c13";"";"c35";"";"c57";"";"c14";"";"c36";"";"c58"]
	xVals = [-4+i for i in 0:8]

	myScatterPlot(title, savePlot, xTitle, yTitle, labels, xVals, c12, xVals, c34, xVals, c56, xVals, c13, xVals, c35, xVals, c57, xVals, c14, xVals, c36, xVals, c58, [-4;-4],[-0.25;0.10])

	J1 = -1
	c12 = MMMfig5helper(J1,J2,1,2)
	c34 = MMMfig5helper(J1,J2,3,4)
	c56 = MMMfig5helper(J1,J2,5,6)
	c13 = MMMfig5helper(J1,J2,1,3)
	c35 = MMMfig5helper(J1,J2,3,5)
	c57 = MMMfig5helper(J1,J2,5,7)
	c14 = MMMfig5helper(J1,J2,1,4)
	c36 = MMMfig5helper(J1,J2,3,6)
	c58 = MMMfig5helper(J1,J2,5,8)

	title = "MMM-fig5-using-MPO-B"

	myScatterPlot(title, savePlot, xTitle, yTitle, labels, xVals, c12, xVals, c34, xVals, c56, xVals, c13, xVals, c35, xVals, c57, xVals, c14, xVals, c36, xVals, c58, [-4;-4],[-0.25;0.10])
end

function MMMfig5helper(Jl, Jr, site1, site2)
	result = [0.0 for i in 1:8]# result = 
	for i in 1:8
		absJ1 = abs(Jl)
		E0, state, sites = dmrgLadder(6, Jl, Jr)
		result[i] = CustOp(sites,state,["Sz";"Sz"],[site1,site2],2)	#Sz, Sx, Sy, all interchangeable as long as both are same
		Jr -= 1
	end
	return result
end

MMMfig5cmp()