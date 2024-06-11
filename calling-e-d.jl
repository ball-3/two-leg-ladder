include("exact-diagonalisation.jl")

let
	for i = 1:12
		print("\n\n For ") 
		print(i) 
		print(" spins:") 
		H = MPO(i)
		#print("\n")
		#display(ED(H))
		print("\nGround State Energy per Spin:\n")
		display(gsE(H)/i)
	end
end
