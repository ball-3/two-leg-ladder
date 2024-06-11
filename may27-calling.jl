using Printf

@printf("helo the file opened\n")

include("heisenberg-one-half-dmrg.jl")
include("ladder-one-half-dmrg.jl")
include("may28-plotting.jl")

include("./semi-inactive/pdiff.jl")
include("./semi-inactive/Hamiltonians.jl")
include("./semi-inactive/exact-diagonalisation.jl")


#TODO this has minimum of two particles to allow for heisenberg 1/2 systems, workaround?
function run(minLength, maxLength, js)
	len = size(js)[1]
	xs = zeros(len)
	ys = zeros(len)
	for j = 1:len
		Jl = js[j,1]
		Jr = js[j,2]
		
		#energies
		energyLDMRG = zeros(maxLength)
		energyLED = zeros(maxLength)
		energyHDMRG = zeros(maxLength)
		energyHED = zeros(maxLength)
		
		#percent differences
		ladderDMRGvED = zeros(maxLength)
		heisenbergDMRGvED = zeros(maxLength)
		DMRGLaddervHeisenberg = zeros(maxLength)
		EDLaddervHeisenberg = zeros(maxLength)
		
		for i = minLength:maxLength
			#find energies
			@printf("This time we are with Jl %d Jr %d and %d spin pairs\n", Jl, Jr, i)
			energyLDMRG[i] = main(i,Jl,Jr)[1]
			energyLED[i] = gsE(Hamiltonians.ladderOneHalf(i,Jl,Jr))
			#energyHDMRG[i] = 2*main2(i,Jl)[1]
			#energyHED[i] = 2*gsE(Hamiltonians.heisenbergOneHalf(i,Jl))	
				
			#find percent differences
			ladderDMRGvED[i] = rdiff(energyLDMRG[i],energyLED[i],i)
			#heisenbergDMRGvED[i] = rdiff(energyHDMRG[i],energyHED[i],i)
			#DMRGLaddervHeisenberg[i] = rdiff(energyLDMRG[i],energyHDMRG[i],i)
			#EDLaddervHeisenberg[i] = rdiff(energyLED[i],energyHED[i],i)
		
			#print energies
			print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Energies ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
			print("Ladder (DMRG): $(energyLDMRG[i])\n")
			print("Ladder   (ED): $(energyLED[i])\n")
			#print("Heise. (DMRG): $(energyHDMRG[i])\n")
			#print("Heisenb. (ED): $(energyHED[i])\n")
			
			#print percent differences
			print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Relative Differences ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
			print("(Ladder)     DMRG vs. ED: $(ladderDMRGvED[i])%\n")
			#print("(Heisenberg) DMRG vs. ED: $(heisenbergDMRGvED[i])%\n")
			#print("(DMRG) Ladder vs. Heise.: $(DMRGLaddervHeisenberg[i])%\n")
			#print("(ED) Ladder vs. Heisenb.: $(EDLaddervHeisenberg[i])%\n")
		end
		
		#plotlin("$j pairs A", true, "x", "y",("Ladder DMRG","Ladder ED","Heis. DMRG","Heis. ED"),energyLDMRG,energyLED,energyHDMRG,energyHED)
		xs[j] = Jl/Jr
		ys[j] = ladderDMRGvED[maxLength]	#this assumes
		
	end
	plot2("graph 2", true, "Jl/Jr", "delta E", ("hmm",),(xs,ys))
end

#m = [10 1;5 1;1 1]

#run(5,5,m)

function myTemp(minPairs,maxPairs,Jl,Jr)
	E = 0
	for i = minPairs:maxPairs
		E = main(i,Jl,Jr)[1]
		diff = rdiff(E/i,-0.75,i)
		@printf("Energy: %f\nDifference: %f\n",E,diff)
	end
end

myTemp(1,50,0,1)
