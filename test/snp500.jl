using Plots, CSV, DataFrames

y = snp500
result = l1tf(y, 50, verbose = true)
x = result.x

plot(ylabel = "log(price)")
plot!(y, label = "S&P500")
plot!(x, label = "l1tf")
png("test/snp500")