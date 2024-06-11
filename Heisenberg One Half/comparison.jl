include("Exact Diagonalisation (mtfishman)/exact_diagonalization.jl")
include("heisenberg-one-half.jl")

function compare(N)
	ed = [NaN, NaN]		#edmrg, vals[1] from exact diagonalisation
	dmrg = NaN

	pdiff1 =NaN
	pdiff2 =NaN

	for j in 1:(N)
		ed = exact_diagonalization.main(N)	#
		dmrg = heisenberg-one-half.main(N)
		pdiff1 = pdiff(ed[1],dmrg)
		pdiff2 = pdiff(ed[2],dmrg)
	
		@printf("For %d sites, exact diagonalisation had E = %f, %f, and DMRG had E = %f.\n
		The differences are %f and %f", N, ed[1], ed[2], dmrg, pdiff1, pdiff2)
	end

	@printf("complete")

	end

function pdiff(expected, actual)
	abs((expected - actual)/expected)*100
end

compare(min(parse(Int64,ARGS[1]),16))	#takes an arg as number of comparisons, if too big, does 16
