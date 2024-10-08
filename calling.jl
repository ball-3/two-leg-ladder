include("./utility-functions/comparison-functions.jl")
include("plotting.jl")

include("./mps/chain.jl")
include("./mps/ladder.jl")
include("./mps/operators.jl")

include("./exact-diagonalisation/Hamiltonians.jl")
include("./exact-diagonalisation/operators.jl")
include("./exact-diagonalisation/exact-diagonalisation.jl")

#displayED(Hamiltonians.ladderOneHalf(2,1,1),5)

function makeRungCouplingPlot(numPairs,gap,JrMax)
	numSites = 2*numPairs
	numDataSets = Int(JrMax/gap)
	title = "average rung pair correlation as a function of coupling strength"
	xTitle = "Jr/Jl"
	yTitle = "rung correlation"
	labels = [""]
	save = true

	Jl = 1
	Jrs = [gap*i for i in 1:numDataSets]
	exp = [0.0 for i in 1:numDataSets]
	
	for i in 1:numDataSets
		energy, state, sites = ladder(numPairs, Jl, Jrs[i])
		for j in 1:2:numPairs
			exp[i] += CustOp(sites,state,["Sz","Sz"],[j,j+1],2)/numSites
		end
	end

	myScatterPlot(title*"$numPairs"*"pairs",save,labels,
	Jrs,exp,
	axisTitles = (xTitle, yTitle)
	)
end

function makeLegCouplingPlot(numPairs,gap,JrMax)
	numSites = 2*numPairs
	numDataSets = Int(JrMax/gap)
	title = "average leg pair correlation as a function of coupling strength"
	xTitle = "Jr/Jl"
	yTitle = "leg correlation"
	labels = [""]
	save = true

	Jl = 1
	Jrs = [gap*i for i in 1:numDataSets]
	exp = [0.0 for i in 1:numDataSets]
	
	for i in 1:numDataSets
		energy, state, sites = ladder(numPairs, Jl, Jrs[i])
		for j in 1:numPairs
			exp[i] += CustOp(sites,state,["Sz","Sz"],[j,j+2],2)/numSites
		end
	end

	myScatterPlot(title*"$numPairs"*"pairs",save,labels,
	Jrs,exp,
	axisTitles = (xTitle, yTitle)
	)
end

function makeMagnetismPlot(numPairs,gap,JrMax)
	numSites = 2*numPairs
	numDataSets = Int(JrMax/gap)
	title = "total magnetism as a function of coupling strength"
	xTitle = "Jr/Jl"
	yTitle = "ground state energy per site"
	labels = [""]
	save = true

	Jl = 1
	Jrs = [gap*i for i in 1:numDataSets]
	exp = [0.0 for i in 1:numDataSets]
	
	for i in 1:numDataSets
		energy, state, sites = ladder(numPairs, Jl, Jrs[i])
		for j in 1:numPairs
			exp[i] += CustOp(sites,state,["Sz"],[j],1)
			#exp += SzOn(state,groundStateWFCT,)
		end
	end

	myScatterPlot(title*"$numPairs"*"pairs",save,labels,
	Jrs,exp,
	axisTitles = (xTitle, yTitle)
	)
end

function makeEnergyPlot(numPairs,gap,JrMax)
	numSites = 2*numPairs
	numDataSets = Int(JrMax/gap)
	title = "ground state energy as a function of coupling strength"
	xTitle = "Jr/Jl"
	yTitle = "ground state energy per site"
	labels = [""]
	save = true

	Jl = 1
	Jrs = [gap*i for i in 1:numDataSets]
	energies = [0.0 for i in 1:numDataSets]

	for i in 1:numDataSets

		energies[i] = ladder(numPairs, Jl, Jrs[i])[1]/(Jrs[i]*numPairs)

	end

	myScatterPlot(title*"$numPairs"*"pairs",save,labels,
	Jrs,energies,
	axisTitles = (xTitle, yTitle)
	)
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

	groundState = ladder(numPairs, Jl, Jr)
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

#TODO brah i need to take an expecation value too otherwise it just sits here ayo
function makeTimeDMRGSzPlot(numPairs, Jl, Jr, t)
	numSites = 2*numPairs
	title = "difference in Magnesitation of state as time evolves vs initial Magnetisation"
	xTitle = "time evolved to"
	yTitle = "Percent Difference (Mz/site)"
	labels = []
	save = true

	DMRGresults = [1.0*i for i in 1:t]
	nTerms = [i for i in 1:t]
	
	for i in 1:t
		psi, sites = timeEvoLadder(2*numPairs,1,1,t,0.1)
		DMRGresults[i] = real(Sz(psi,sites))
	end

	pDiffs = [pdiff(DMRGresults[i],DMRGresults[1]) for i in 1:t]

	myScatterPlot(title*" $numPairs"*" pairs",save,labels,
	nTerms, pDiffs,
	axisTitles = (xTitle, yTitle))
end


makeTimeDMRGSzPlot(5,1,1,25)

#makeCorrelationPlot(2,1,10)