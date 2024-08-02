include("plotting.jl")

x = [-3;-2;-1;0;1;2;3]
y = [1;2;3;4;5;6;7]

myScatterPlot("temp-testing->:0", true, "x", "y", ["series1","","series2","","series3"], x, y, x, 2*y,[-3;-3],[-16;0])