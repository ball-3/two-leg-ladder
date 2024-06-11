using Plots

x = range(0,10, length=100)
y = sin.(x)
z = #data
plot(x,y)
plot(x, z,
