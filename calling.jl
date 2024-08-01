include("./utility-functions/comparison-functions.jl")
include("plotting.jl")

include("./dmrg/chain.jl")
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

function callingABCTest(nPairs, doED::Bool, doMPS::Bool, Jl, Jr, site1, site2)
	sigma = 0.5*[σx,σy,σz]
	sigmaStr = ["Sx","Sy","Sz"]
	
	expA, expB, expC = 0, 0, 0

	if(doED)
		groundState = ED(Hamiltonians.ladderOneHalf(nPairs,Jl,Jr))[2][1]

		for j in 1:3
			op = [sigma[j],sigma[j]]
			opStr = "<"*sigmaStr[j]*sigmaStr[j]*">"

			expA = CustOp(groundState,op,[site1,site2])
			#expB = CustOp(groundState,op,[i,i+2])
			#expC = CustOp(groundState,op,[i+1,i+3])

			exp = expA - 0.5*expB - 0.5*expC
			#@printf("Using ED for operator %s  : %f = %f - 1/2 %f - 1/2 %f\n", opStr, exp, expA, expB, expC)

		end
	end

	if(doMPS)
		groundState = dmrgLadder(nPairs, Jl, Jr)
		groundStateWFCT = groundState[2]
		sites = groundState[3]

		for j in 1:3
			op = [sigmaStr[j],sigmaStr[j]]
			opStr = "<"*sigmaStr[j]*sigmaStr[j]*">"

			expAMPO = CustOp(sites,groundStateWFCT,op,[site1,site2],2)
			return [expA, expAMPO]

			if (nPairs > 1)
				#expB = CustOp(sites,groundStateWFCT,op,[i,i+2],2)
				#expC = CustOp(sites,groundStateWFCT,op,[i+1,i+3],2)
			end

			#TODO this is a temp change to fit my plot fct, edit later?
			return [expA,expB,expC]
		end
	end
end

function makeABCPlot(numPairs,Jl,Jr)
	numDataSets = 16
	title = "correlations by site pair"
	xTitle = "site pair"
	yTitle = "spin-spin correlation"
	save = true
	
	labels = ["c12-e";"";"c12-m";"";"c21-e";"";"c21-m";"";"c34-e";"";"c34-m";"";"c43-e";"";"c43-m";"";"c56-e";"";"c56-m";"";"c65-e";"";"c65-m";"";"c78-e";"";"c78-m";"";"c87-e";"";"c87-m";"";]
	c12e, c12m = callingABCTest(numPairs,true,true,Jl,Jr,1,2)
	c21e, c21m = callingABCTest(numPairs,true,true,Jl,Jr,2,1)
	c34e, c34m = callingABCTest(numPairs,true,true,Jl,Jr,3,4)
	c43e, c43m = callingABCTest(numPairs,true,true,Jl,Jr,4,3)
	c56e, c56m = callingABCTest(numPairs,true,true,Jl,Jr,5,6)
	c65e, c65m = callingABCTest(numPairs,true,true,Jl,Jr,6,5)
	c78e, c78m = callingABCTest(numPairs,true,true,Jl,Jr,7,8)
	c87e, c87m = callingABCTest(numPairs,true,true,Jl,Jr,8,7)

	myScatterPlot(title*"A",save,xTitle,yTitle,labels,
	[1], [Float64(c12e)],
	[1], [Float64(c12m)],
	[1], [Float64(c21e)],
	[1], [Float64(c21m)],
	[2], [Float64(c34e)],
	[2], [Float64(c34m)],
	[2], [Float64(c43e)],
	[2], [Float64(c43m)],
	[3], [Float64(c56e)],
	[3], [Float64(c56m)],
	[3], [Float64(c65e)],
	[3], [Float64(c65m)],
	[4], [Float64(c78e)],
	[4], [Float64(c78m)],
	[4], [Float64(c87e)],
	[4], [Float64(c87m)]
	)
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

function terrible(val::Bool)
	return val
end

#callingABCTest(1,terrible(false),terrible(false),1,1)
#callingABCTest(1,false,false,1,1)

#=
groundState = ED(Hamiltonians.ladderOneHalf(4,1,1))[2][1]
c12 = CustOp(groundState,[σz,σz],[8,1])
c34 = CustOp(groundState,[σz,σz],[4,3])
c56 = CustOp(groundState,[σz,σz],[6,5])
c78 = CustOp(groundState,[σz,σz],[8,7])
c58 = CustOp(groundState,[σz,σz],[5,8])
c67 = CustOp(groundState,[σz,σz],[6,7])

print("C12: $(c12)  C34: $(c34)  C56: $(c56)  C78: $(c78)  C58: $(c58)  C67: $(c67)\n")
=#
#=
groundState = dmrgLadder(4, 1, 1)
		groundStateWFCT = groundState[2]
		sites = groundState[3]
c12 = CustOp(sites,groundStateWFCT,["Sz","Sz"],[1,2],2)
c34 = CustOp(sites,groundStateWFCT,["Sz","Sz"],[3,4],2)
c56 = CustOp(sites,groundStateWFCT,["Sz","Sz"],[5,6],2)
c78 = CustOp(sites,groundStateWFCT,["Sz","Sz"],[7,8],2)
c58 = CustOp(sites,groundStateWFCT,["Sz","Sz"],[5,5],2)
c67 = CustOp(sites,groundStateWFCT,["Sz","Sz"],[6,7],2)

print("C12: $(c12)  C34: $(c34)  C56: $(c56)  C78: $(c78)  C58: $(c58)  C67: $(c67)\n")
=#

makeABCPlot(4,1,1)