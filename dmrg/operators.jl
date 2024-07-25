function CustOp(sites, state, operators, sitesApplied, prodNum)	
	numSitesApplied = length(sitesApplied)
	numOperators = length(operators)
	
	if numOperators > numSitesApplied
		print("WARNING: one or more operators specified have no state to apply to and will be ignored\n")
	end
	
	print("This function is for testing and doesnt actually work yet\n")
	os = OpSum()
	for j in 1:(numSitesApplied-1)

		#for k in 1:prodNum
			#inter dimer
			os += operators[j], sitesApplied[j]
			os *= operators[j+1], sitesApplied[j+1]
		#end
	end

	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal
end

function S2(sites, state)	
	numSites = length(sites)
	
	os = OpSum()
	for j in 1:numSites
	#os += "Sx",1,"Sx",1
	os += "Sx * Sx",j
	#os += "Sy",1,"Sy",1
	os += "Sy * Sy",j
	#os += "Sz",1,"Sz",1
	os += "Sz * Sz",j
	end
	
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)/length(sites)
	
	return expVal
end

function Sz(sites, state)	
	numSites = length(sites)
	
	os = OpSum()
	for j in 1:(numSites)
	os += "Sz",j
	end
	
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal/length(sites)
end

function Sx(sites, state)	
	numSites = length(sites)
	
	os = OpSum()
	for j in 1:(numSites)
	os += "Sx",j
	end
	
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal/length(sites)
end

function Sy(sites, state)	
	numSites = length(sites)
	
	os = OpSum()
	for j in 1:(numSites)
	os += "Sy",j
	end
	
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal/length(sites)
end

function Splus(sites, state)	
	numSites = length(sites)
	
	os = OpSum()
	for j in 1:(numSites)
	os += "S+",j
	end
	
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal/length(sites)
end

function Sminus(sites, state)	
	numSites = length(sites)
	
	os = OpSum()
	for j in 1:(numSites)
	os += "S-",j
	end
	
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal/length(sites)
end
