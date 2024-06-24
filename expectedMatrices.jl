using Plots

lanczos = [ NaN 0.0 0.2 0.4 0.6 0.8 1.0;
	   4 -0.500000 -0.503782 -0.515856 -0.536654 -0.565808 -0.602511;
	   6 -0.467133 -0.472005 -0.487058 -0.511605 -0.544412 -0.584437;
	   8 -0.456387 -0.462023 -0.478718 -0.504825 -0.538967 -0.580203;
	   10 -0.451545 -0.457740 -0.475373 -0.502232 -0.537031 -0.578860;
	   12 -0.448949 -0.455583 -0.473761 -0.501045 -0.536229 -0.578375;
	   24 -0.444584 NaN NaN NaN NaN NaN;
	   32 NaN NaN NaN NaN NaN NaN;
	   Inf -0.443147 NaN NaN NaN NaN NaN ]
	   
montecarlo = [  NaN 0.0 0.2 0.4 0.6 0.8 1.0;
	   4 NaN NaN NaN NaN NaN NaN;
	   6 NaN NaN NaN NaN NaN NaN;
	   8 -0.4563 -0.4623 -0.4790 -0.5047 -0.5390 -0.5800;
	   10 NaN NaN NaN NaN NaN NaN;
	   12 NaN NaN NaN NaN NaN NaN;
	   24 -0.4447 -0.4526 -0.4719 -0.4995 -0.5354 -0.5781;
	   32 NaN NaN -0.4726 -0.5001 -0.5357 -0.5784;
	   Inf -0.443 -0.45 -0.472 -0.500 -0.535 -0.578 ]
	   
include("ladder-one-half-dmrg.jl")
include("./semi-inactive/pdiff.jl")
	   
function compare()
	vsLanczos = zeros(9,7)
	vsMC = zeros(9,7)
	
	#setting up column and row labels
	for i = 1:7
		vsLanczos[1,i] = lanczos[1,i]
		vsMC[1,i] = montecarlo[1,i]
	end
	for i = 2:9
		vsLanczos[i,1] = lanczos[i,1]
		vsMC[i,1] = montecarlo[i,1]
	end

	for i = 2:4#should be nine but ignoring infinite case
		for j = 2:7
		#i is num of spin pairs
		#E0/2LJ
			numPairs = lanczos[i,1]
			ratio = lanczos[1,j]
			dmrg = dmrgRES(numPairs,ratio)
			print("$(i) and $(j) , ratio; $(ratio), numPairs = $(numPairs)\n")
			#print("Lanczos: $(lanczos[i,j]) DMRG: $(dmrg)\n")
			#vsLanczos[i,j] = rdiff(lanczos[i,j],dmrg,numPairs*2)
			#vsMC[i,j] = rdiff(montecarlo[i,j],dmrg,numPairs*2)
			EperSite = -(3/8)*(ratio)-(3/16)*(1/ratio)
			vsLanczos[i,j] = rdiff(dmrg,EperSite,numPairs*2)
			vsMC[i,j] = rdiff(dmrg,EperSite,numPairs*2)
		end
	end
	print("Lanczos:\n")
	display(vsLanczos)
	print("Monte Carlo:\n")
	display(vsMC)

end

function plotComparison(col)
	x = [4;6;8;10;12;24]
	yLanc = [lanczos[i,col] for i in 2:8]
	yMC = [montecarlo[i,col] for i in 2:8]
	ratio = lanczos[1,col]
	val = -(3/8)*(ratio)-(3/16)*(10000/ratio)
	yCalc = [val;val;val;val;val;val;val]
	
	yDMRG = [dmrgRES(i,ratio) for i in 2:8]

	plot()	
	scatter!(x,yLanc, ms = 10, label="lanczos results")
	scatter!(x,yMC, ms = 7, label="montecarlo results")
	scatter!(x,yDMRG, ms = 5, label="dmrg results")
	if col != 2
		scatter!(x,yCalc,label="function results")
	end
	
	png("results for jr = $(ratio) and jl = 1 attempt two")
end

function dmrgRES(pairs, ratio)
	jl = 1
	jr = jl*ratio
	return (dmrgLadder(pairs, jl, jr,-1)[1])/(2*pairs*jl)
end

for i = 2:7
	plotComparison(i)
end
