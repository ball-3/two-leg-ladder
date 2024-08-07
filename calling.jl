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

#designed for an odd number of spin pairs
function makeCorrelationPlot(numPairs,Jl,Jr)
	numSites = 2*numPairs
	title = "correlations by site pair"
	xTitle = "site pair"
	yTitle = "spin-spin correlation"
	save = true

	#correlations
	#	c	i  i+2 i+4 i+6  .....		c	i  i+2 i+4 i+6  .....	c	i   i+2  i+4  i+6  .....	c	i   i+2  i+4  i+6  .....	
	#	D	*	*	*	*				U	*	*	*	*			R	* -> * -> * -> *			L	*	 *	  *	   *			
	#	o	|	|	|	|	.....		p	^	^	^	^			i								e
	#	w	v	v	v	v				 	|	|	|	|	.....	g					   .....	f					   .....	
	#	n	*	*	*	*					*	*	*	*			ht	*	 *	  *	   *			t	* <- * <- * <- *			
	cDown = [0.0 for i in 1:numPairs]
	cUp = [0.0 for i in 1:numPairs]
	cLeft = [0.0 for i in 1:(numPairs-1)]
	cRight = [0.0 for i in 1:(numPairs-1)]

	#y\ =\ -\operatorname{abs}\left(\frac{l}{2}-x+\frac{1}{2}\right)+\ \frac{l}{2}+\frac{1}{2}
	hlen = numPairs/2
	if iseven(numPairs)
		hlen += 1/2
	end
	xVals1 = [(-abs(hlen - i) + hlen) for i in 1:numPairs]
	xVals2 = [(-abs(hlen - i) + hlen) for i in 1:(numPairs-1)]

	groundState = dmrgLadder(numPairs, Jl, Jr)
	groundStateWFCT = groundState[2]
	sites = groundState[3]
	op = ["Sz","Sz"]
	
	labels = ["correlation down";"correlation up";"correlation left top";"correlation right bottom"]
	for i in 1:numPairs
		cDown[i] = CustOp(sites,groundStateWFCT,op,[i,i+1],2)
		cUp[i] = CustOp(sites,groundStateWFCT,op,[i+1,i],2)
	end
	ctr = 0
	for i in 1:2:(numSites-3)
		ctr += 1
		cRight[ctr] = CustOp(sites,groundStateWFCT,op,[i,i+2],2)
		cLeft[ctr] = CustOp(sites,groundStateWFCT,op,[i+3,i+1],2)
	end

	myScatterPlot(title*"$numPairs"*"pairs",save,labels,
	xVals1, cDown,
	xVals1, cUp,
	xVals2, cLeft,
	xVals2, cRight,
	axisTitles = (xTitle, yTitle)
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

	myScatterPlot(title, savePlot, labels, xVals, yVals1, xVals, yVals2, axisTitles = (xTitle, yTitle), axisLims = ([-4;4],[-18;2]))
end

function MMMfig4helper(Jl, Jr)
	result = [0.0 for i in 1:7]# result = ground state energy by absolute value of Jl, y values in fig 4
	for i in 1:7
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
	labels = ["c12";"c34";"c56";"c13";"c35";"c57";"c14";"c36";"c58"]
	xVals = [-4+i for i in 0:8]

	myScatterPlot(title, savePlot, labels, xVals, c12, xVals, c34, xVals, c56, xVals, c13, xVals, c35, xVals, c57, xVals, c14, xVals, c36, xVals, c58, axisTitles = (xTitle, yTitle), axisLims = ([-5;5],[-0.30;0.15]))

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

	myScatterPlot(title, savePlot, labels, xVals, c12, xVals, c34, xVals, c56, xVals, c13, xVals, c35, xVals, c57, xVals, c14, xVals, c36, xVals, c58, axisTitles = (xTitle, yTitle), axisLims = ([-5;5],[-0.30;0.15]))

end

function MMMfig5helper(Jl, Jr, site1, site2)
	result = [0.0 for i in 1:9]# result = 
	for i in 1:9
		absJ1 = abs(Jl)
		E0, state, sites = dmrgLadder(6, Jl, Jr)
		result[i] = CustOp(sites,state,["Sz";"Sz"],[site1,site2],2)	#Sz, Sx, Sy, all interchangeable as long as both are same
		Jr -= 1
	end
	return result
end