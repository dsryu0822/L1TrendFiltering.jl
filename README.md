# L1TrendFiltering.jl

A julia implementation for [l1 trend filtering](https://web.stanford.edu/~gorin/papers/l1_trend_filter.pdf). It's just a translation of [original matlab code](https://web.stanford.edu/~boyd/l1_tf/), but something minor has been changed and readability is improved.

## Usage

```julia
using L1TrendFiltering
using Plots, CSV, DataFrames

y = snp500
result = l1tf(y, 50, verbose = true)
x = result.x

plot(ylabel = "log(price)")
plot!(y, label = "S&P500")
plot!(x, label = "l1tf")
```

![snp](test/snp500.png)
