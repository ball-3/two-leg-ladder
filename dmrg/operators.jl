function CustOp(sites, state, operators, sitesApplied, prodNum)	
	numSitesApplied = length(sitesApplied)
	numOperators = length(operators)
	
	if numOperators > numSitesApplied
		print("WARNING: one or more operators specified have no state to apply to and will be ignored\n")
	end
	
	os = OpSum()
	#TODO the following needs thorough testing
	#=
	for j in 1:prodNum:(numSitesApplied)
		
		os += operators[j], sitesApplied[j]

		for k in 1:(prodNum-1)
			os *= operators[j+k], sitesApplied[j+k]
		end
	end
	=#
	if prodNum == 2
		for j in 1:Integer(numSitesApplied/2)
			os += operators[j], sitesApplied[j]
			os *= operators[j+1], sitesApplied[j+1]
		end
	else#TODO this assumes prodnum == 1 'B)
		for j in 1:Integer(numSitesApplied)
			os += operators[j], sitesApplied[j]
		end
	end
	op = MPO(os, sites)
	
	expVal = inner(state',op,state)
	
	return expVal
end

#finds the expectation value of S^2 operator given a state
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

#finds the expectation value of S_z operator given a state
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

#finds the expectation value of S_x operator given a state
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

#finds the expectation value of S_y operator given a state
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

#finds the expectation value of S^+ operator given a state
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

#finds the expectation value of S^- operator given a state
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