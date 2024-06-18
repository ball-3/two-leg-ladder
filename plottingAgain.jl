using Plots

include("ladder-one-half-dmrg.jl")

function plot3(nVals,jlVals,jrVals)
	
	cols = length(nVals)
	rows = length(jlVals)
	jrs = length(jrVals)
	p = Matrix{Plots.Plot{Plots.GRBackend}}(undef, rows, cols)
	
	for jl = 1:rows
		for n =1:cols
			E = [dmrgLadder(n,jl,i)[1] for i in 1:jrs]
			p[jl,n] = plot3sub([jl,[jrVals[1];jrVals[jrs]],n], jrVals, E)
		end
	end			
	plot(p[1,1],p[1,2],p[1,3],p[1,4],p[1,5],p[2,1],p[2,2],p[2,3],p[2,4],p[2,5],p[3,1],p[3,2],p[3,3],p[3,4],p[3,5],p[4,1],p[4,2],p[4,3],p[4,4],p[4,5],p[5,1],p[5,2],p[5,3],p[5,4],p[5,5],p[6,1],p[6,2],p[6,3],p[6,4],p[6,5],p[7,1],p[7,2],p[7,3],p[7,4],p[7,5], layout=(rows,cols), legend = false)
	
	png("hewwop")
end

function plot3sub(fData,xs,ys)
	p = plotFunction(expectedEnergy,fData[1],fData[2],fData[3])
	scatter!(p, xs,ys)
	
	return p
end

#expect bounds like [[xmin xmax];[ymin ymax]]
function plotFunction(f,Jl, bounds, N, stepGap=100)
	x = [i for i in bounds[1]:stepGap:bounds[2]]
	y = f.(x,Jl,N) 
	return plot!(x,y)
end

function expectedEnergy(Jr,Jl,N)
	Egs = (-0.75*Jr-0.375*Jl*Jl/Jr)*N
end


nVals = [1;10;100;1000;10000]
jrVals = [0.1;1;10;100;1000;10000]
jlVals = [0.0001;0.001;0.01;0.1;1;10;100]

plot3(nVals,jlVals,jrVals)
