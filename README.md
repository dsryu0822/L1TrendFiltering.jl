# L1TrendFiltering.jl

A julia implementation for [l1 trend filtering](https://web.stanford.edu/~gorin/papers/l1_trend_filter.pdf). It's just a translation of [original matlab code](https://web.stanford.edu/~boyd/l1_tf/), but something minor has been changed and readability is improved.

## Usage

 - This package includes a dataset `snp500` as a vector.
 - The function `l1tf` returns the result `L1tf` object. `L1tf` has two properties:
   - `x`: `Vector{T}`, output of $l_{1}$ trend filtering.
   - `solved`: `Bool`, if the solution is converged, then `true`.

In this package, `l1tf` allows changing other parameters, for instance, linesearch parameters $\alpha$ and $\beta$.

### Example

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
