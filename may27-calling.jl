include("heisenberg-one-half-dmrg.jl")
include("ladder-one-half-dmrg.jl")
include("may28-plotting.jl")
include("expectedMatrices.jl")

include("./semi-inactive/pdiff.jl")
include("./semi-inactive/Hamiltonians.jl")
include("./semi-inactive/exact-diagonalisation.jl")


#TODO this has minimum of two particles to allow for heisenberg 1/2 systems, workaround?
#note that this function has outdated method names and will not compile

#dmrgLadder(numPairs, jl, jr, operator=-1)

function run(minPairs, maxPairs, js)
	len = size(js)[1]
	xs = zeros(len)
	ys = zeros(len)
	for j = 1:len
		Jl = js[j,1]
		Jr = js[j,2]
		
		#energies
		energyLDMRG = zeros(maxPairs)
		energyLED = zeros(maxPairs)
		energyHDMRG = zeros(maxPairs)
		energyHED = zeros(maxPairs)
		
		#percent differences
		ladderDMRGvED = zeros(maxPairs)
		heisenbergDMRGvED = zeros(maxPairs)
		DMRGLaddervHeisenberg = zeros(maxPairs)
		EDLaddervHeisenberg = zeros(maxPairs)
		ladderLANCvED = zeros(maxPairs)
		ladderLANCvDMRG = zeros(maxPairs)
		
		for i = minPairs:2:maxPairs
			#find energies
			@printf("This time we are with Jl %f Jr %f and %d spin pairs\n", Jl, Jr, i)
			energyLDMRG[i] = dmrgLadder(i,Jl,Jr,-1)[1]/(2*i*Jl)
			energyLED[i] = gsE(Hamiltonians.ladderOneHalf(i,Jl,Jr))/(2*i*Jl)
			#energyHDMRG[i] = 2*main2(i,Jl)[1]
			#energyHED[i] = 2*gsE(Hamiltonians.heisenbergOneHalf(i,Jl))	
				
			#find percent differences
			# save the turtles :D
			ladderDMRGvED[i] = rdiff(energyLDMRG[i],energyLED[i],i)
			ladderLANCvDMRG[i] = rdiff(lanczos[i+1,j+1],energyLDMRG[i],i)
			ladderLANCvED[i] = rdiff(lanczos[i+1,j+1],energyLED[i],i)
			#heisenbergDMRGvED[i] = rdiff(energyHDMRG[i],energyHED[i],i)
			#DMRGLaddervHeisenberg[i] = rdiff(energyLDMRG[i],energyHDMRG[i],i)
			#EDLaddervHeisenberg[i] = rdiff(energyLED[i],energyHED[i],i)
		
			#print energies
			@printf("This time we are with Jl %f Jr %f and %d spin pairs\n", Jl, Jr, i)
			print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Energies ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
			print("Ladder (DMRG): $(energyLDMRG[i])\n")
			print("Ladder   (ED): $(energyLED[i])\n")
			print("Ladder (Lanc): $(lanczos[i+1,j+1])\n")
			#print("Heise. (DMRG): $(energyHDMRG[i])\n")
			#print("Heisenb. (ED): $(energyHED[i])\n")
			
			#print percent differences
			print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Relative Differences ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
			print("(Ladder)     DMRG vs. ED: $(ladderDMRGvED[i])%\n")
			print("(Ladder)   LANC vs. DMRG: $(ladderLANCvDMRG[i])%\n")
			print("(Ladder)     LANC vs. ED: $(ladderLANCvED[i])%\n")
			#print("(Heisenberg) DMRG vs. ED: $(heisenbergDMRGvED[i])%\n")
			#print("(DMRG) Ladder vs. Heise.: $(DMRGLaddervHeisenberg[i])%\n")
			#print("(ED) Ladder vs. Heisenb.: $(EDLaddervHeisenberg[i])%\n")
		end
		
		#plotlin("$j pairs A", true, "x", "y",("Ladder DMRG","Ladder ED","Heis. DMRG","Heis. ED"),energyLDMRG,energyLED,energyHDMRG,energyHED)

	end
end

function cmpEgsN(minPairs,maxPairs,jVals)
	len = size(jVals)[1]
	JrVals = zeros(len)
	JlVals = zeros(len)
	EDgse = zeros(len)
	DMRGgse = zeros(len)
	EXPgse = zeros(len)
	for i = 1:len
		JrVals[i] = jVals[i,2]
		JlVals[i] = jVals[i,1]
	end
	for i = minPairs:maxPairs	
		for j = 1:len
			Jl = JlVals[j]
			Jr = JrVals[j]
			#EDgse[j] = gsE(Hamiltonians.ladderOneHalf(i,Jl,Jr))/(2*i)
			try
				DMRGgse[j] = dmrgLadder(i,Jl,Jr)[1]/(2*i)
				EXPgse[j] = (-0.75*Jr-0.375*Jl*Jl/Jr)/2
			catch e
				DMRGgse[j] = 0
				EXPgse[j] = 0
			end
		end
		#plot2("Ground state energy versus Jr strength for $(i) spin pairs, large range",true,"Jr","Egs",[],i,JrVals,JlVals,Egs)
		plot2("Ground state energy versus Jr strength for $(i) spin pairs, ED and DMRG",true,"Jl","Egs",["Function";"Jl";"DMRG"],i,JlVals,EXPgse,JlVals,DMRGgse)
		#plot2(title::String, save::Bool, xAxisTitle, yAxisTitle, labels, numPairs, Jl, data...)
	end
end

jr = 1
jlVals = [jl for jl in 0.01:0.01:0.1]
jVals = zeros(length(jlVals),2)

for i in 1:length(jlVals)
	jVals[i,1] = jlVals[i]
	jVals[i,2] = jr
end

cmpEgsN(1,10,jVals)
#run(4,4,[1 0;1 0.2;1 0.4;1 0.6;1 0.8;1 1])
