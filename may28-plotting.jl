#Note : determined errors when x, yn are not all of the same type (row OR column vector), so must ensure same type
#ensure labels[i] are strings to refer to the ith 
using Plots

#TODO make this function more modular

#defaultcolours = [red parse(Colorant, "blue1") green parse(Colorant, "blue1") blue parse(Colorant, "blue1")]

function plotlin(title::String, save::Bool, xAxisTitle::String, yAxisTitle::String,labels=(),#=colours=defaultcolours,=# y...)

	if (isempty(y)) return end
	
	ys = size(y,1)
	
	if ys > size(labels,1)
		ls = size(labels,1)
		newlabels = String[]
		#filling empty labels
		for i = 1:ls
			push!(newlabels,labels[i])
		end
		for i = ls:ys
			push!(newlabels, "untitled$(i+1)")	#this probably shouldnt have a +1, where is the logic error
		end
		labels = newlabels
	end
	#colours = defaultcolours
	#colourslength =  max(size(colours)[1],size(colours)[2])
	
	plot() #TODO what is supposed to be here?
	plot!(legend=:topright)
	title!(title)
	xlabel!(xAxisTitle)
	ylabel!(yAxisTitle)
	
	for i = 1:ys
		n = size(y[i],1)
		x = [1/(2*j) for j in 1:n]
		scatter!(x,y[i],label=labels[i])
		#scatter!(x, y[i], label="y$i", mc=:colours[(i%colourslength)], ms=(ys-i), ma=0.5)
	end

	if save
		png(title)
	end
end

function plotlog()
		#plot!(xscale=:log10, yscale=:log10, minorgrid=true)
end

function plot2(title::String, save::Bool, xAxisTitle, yAxisTitle, labels, numPairs, data...)
	
	if (isempty(data)) return end
	
	ndata = size(data,1)
	
	if ndata > size(labels,1)
		ls = size(labels,1)
		newlabels = String[]
		#filling empty labels
		for i = 1:ls
			push!(newlabels,labels[i])
		end
		for i = ls:ndata
			push!(newlabels, "untitled$(i+1)")	#this probably shouldnt have a +1, where is the logic error
		end
		labels = newlabels
	end
	
	plot()
	plot!(legend=:topright)
	title!(title)
	xlabel!(xAxisTitle)
	ylabel!(yAxisTitle)
	
	xmin = 0
	xmax = 0
	
	for i = 1:2:ndata
		xVals = data[i]
		yVals = data[i+1]	#*
		
		scatter!(xVals, yVals,label=labels[i])
		for j = 1:(min(length(xVals),length(yVals)))
			x = xVals[j]
			y = yVals[j]
		
			xmin = min(x,xmin)
			xmax = max(x,xmax)
		end
		plotFunction(expectedEnergy,[xmin;xmax],numPairs)
		if save
			png(title)
		end
	
	end	
end

#expect bounds like [[xmin xmax];[ymin ymax]]
function plotFunction(f, bounds, N, stepGap=100)
	x = [i for i in bounds[1]:stepGap:bounds[2]]
	y = f.(x,N) 
	plot!(x,y)
end

function expectedEnergy(Jr,N)
	Egs = -0.75*Jr*N
end
