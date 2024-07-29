using Plots

markershapes = [:circle;:circle;:star5;:star5;:rect;:rect;:h;:h;:star7;:star7;:ltriangle;:ltriangle;:star4;:star4;:utriangle;:utriangle;:dtriangle;:dtriangle]

function myScatterPlot(title::String, save::Bool, xAxisTitle, yAxisTitle, labels, data...)
	
	if (isempty(data)) return end
	
	#number of data series
	nSeries = size(data,1)
	
	if nSeries > size(labels,1)
		ls = size(labels,1)
		newlabels = String[]
		#filling empty labels
		for i = 1:ls
			push!(newlabels,labels[i])
		end
		for i = ls:nSeries
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
	
	for i = 1:2:nSeries
		xVals = data[i]
		yVals = data[i+1]
		
		if (nSeries < 6)
			scatter!(xVals, yVals,ms = 2*nSeries-2*i,label=labels[i])
		else
			scatter!(xVals, yVals,ms = 4, shape = markershapes[i%(length(markershapes))] ,label=labels[i])
		end

		for j = 1:(min(length(xVals),length(yVals)))
			x = xVals[j]
			y = yVals[j]
		
			xmin = min(x,xmin)
			xmax = max(x,xmax)
		end
		#plotFunction(expectedEnergy,[xmin;xmax],numPairs)
		if save
			png(title)
		end
	
	end	
end

#expect bounds like [[xmin xmax];[ymin ymax]]
function plotFunction(f, bounds, N, stepGap=100)
	jl = [i for i in bounds[1]:stepGap:bounds[2]]
	y = f.(1,jl,N) 
	plot!(jl,y)
end

function expectedEnergy(Jr,Jl,N)
	Egs = (-0.75*Jr-0.375*Jl*Jl/Jr)/2
end



