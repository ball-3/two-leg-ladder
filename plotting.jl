using Plots

markershapes = [:circle;:circle;:star5;:star5;:rect;:rect;:h;:h;:star7;:star7;:ltriangle;:ltriangle;:star4;:star4;:utriangle;:utriangle;:dtriangle;:dtriangle]

function myScatterPlot(title::String, save::Bool, labels, data...; axisTitles = "auto", axisLims = "auto", axisTicks = "auto")
	
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
	
	plot(dpi = 600)
	#plot!(legend=:topright)
	title!(title)
	
	if !(typeof(axisTitles) == String && axisTitles == "auto")
		plot!(xlabel = axisTitles[1], ylabel = axisTitles[2])
	end
	if !(typeof(axisLims) == String && axisLims == "auto")
		plot!(xlims = axisLims[1], ylims = axisLims[2])
	end
	if !(typeof(axisTicks) == String && axisTicks == "auto")
		plot!(xticks = axisTicks[1], yticks = axisTicks[2])
	end
	
	xmin = 0
	xmax = 0

	counter = 1
	
	for i = 1:2:(nSeries-1)
		xVals = data[i]
		yVals = data[i+1]

		plot!(xVals, yVals,ms = 6, shape = markershapes[i%(length(markershapes))] ,label=labels[counter], alpha = 0.5)

		counter += 1

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



