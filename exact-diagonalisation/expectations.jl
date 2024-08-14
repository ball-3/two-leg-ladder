include("kronecker-product.jl")
include("pauli-matrices.jl")
include("exact-diagonalisation.jl")	#do i need this?

#get state by exact diagonalisation
#display(ED(Hamiltonians.ladderOneHalf(2,1,1))[2][1])
#so each function takes in ED hamiltonian hamtype (params) [2] [1] (this is column vector/ ket of ground state)
#and creates a bra (ket')

function normalisationTest(state)
	ket = state
	bra = state'
	return bra*ket
end

function CustOp(state, operators, sitesApplied)
	len = length(state)
	dim = Int(log2(len))
	ket = state
	bra = state'
	
	op = buildOperator(operators,sitesApplied,dim)
	
	return (bra*op*ket)
end

function Sz(state)
	len = length(state)
	dim = Int(log2(len))
	ket = state
	bra = state'
	
	op = buildOperator([σz], "all", dim)
	
	return (bra*op*ket)
end

function SzOn(state, siteApplied)
	len = length(state)
	dim = Int(log2(len))
	ket = state
	bra = state'
	siteApplied = [siteApplied]
	
	op = buildOperator([σz], siteApplied, dim)
	
	return (bra*op*ket)
end

function SyOn(state, siteApplied)
	len = length(state)
	dim = Int(log2(len))
	ket = state
	bra = state'
	siteApplied = [siteApplied]
	
	op = buildOperator([σy], siteApplied, dim)
	
	return (bra*op*ket)
end

#TODO do i need to multiply by one half as im using pauli matrices or should i incorporate that in "operators" matrix
function buildOperator(operators, sitesApplied, dim)
	
	Op = 1
	
	if ((typeof(sitesApplied) == String) && (lowercase(sitesApplied) == "all"))
		# ~~~~ applying operators to all sites, unspecified operators as identities ~~~~
		len = min(length(operators), dim)
		for i in 1:len
			Op = kron(Op, operators[i])
		end
		for i in (len+1):dim
			Op = kron(Op, I2)
		end
		return Op
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	else
		lenSites = length(sitesApplied)
		# **** checking that all sitesApplied have an operator, and if not padding with identities ****
		lenOps = min(length(operators), dim)
		if lenOps < lenSites
			operators2 = [I2 for i in 1:dim]
			for i in 1:lenOps
				operators2 = operators[i]
			end
			for i in (lenOps+1):dim
				operators2 = I2
			end
		end
		# ****************************************************************************************************
		counter = 1
		# $$$$ applying operators to sitesApplied $$$$
		for i in 1:dim
				if ((counter > lenSites) || (i < sitesApplied[counter]))
					Op = kron(Op, I2)
				else
					Op = kron(Op,operators[counter])
					counter = counter+1
				end
		end
		return Op
		# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	end
end

#currently for only ground state
function timeEvolution(tMax, tGap, H)
	diag = ED(H)
	E = diag[1][1]
	state0 = diag[2][1]
	dim = length(state0)
	state = [0.0 + 0.0*im for i in 1:dim]
	state += state0

	for t in 0:tGap:tMax
		x = -t*E
		state *= (cos(x) + im*sin(x))
		state = normalize(state)
		#normalise the state here and now
	end
	
	return state
end

H0 = [1 0 0 0; 0 -1 1 0; 0 1 -1 0; 0 0 0 1]
c = timeEvolution(5,0.1, H0)
display(c)
display(normalisationTest(c))