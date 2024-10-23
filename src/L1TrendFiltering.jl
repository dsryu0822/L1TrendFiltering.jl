module L1TrendFiltering

    using LinearAlgebra, SparseArrays, Printf

    struct L1tf
        x::Vector{Float64}
        solved::Bool
    end
    function Base.show(io::IO, result::L1tf)
        print(io, "$(length(result.x))-element l1tf result, solved: $(result.solved)")
    end

    """
        l1tf(y, λ; α = 1e-2, β = 5e-1, μ = 2, MAXITER = 40, MAXLSITER = 20, atol = 1e-4, verbose = false)

        Find the solution of the L1 trend filtering problem using interior-point method:
        minimize ``{\\frac{1}{2}} \\left\\| y - x \\right\\|_{2}^{2} + \\lambda \\left\\| D x \\right\\|_{1}``
        where ``D`` is the second difference matrix.

        The implmentation is based on original MATLAB code by authors.
        
        ### Optional Arguments
        - `α = 1e-2`: backtracking linesearch parameter (0,0.5]
        - `β = 5e-1`: backtracking linesearch parameter (0,1)
        - `μ = 2`: IPM parameter: t update
        - `MAXITER = 40`: IPM parameter: max iteration of IPM
        - `MAXLSITER = 20`: IPM parameter: max iteration of line search
        - `atol = 1e-4`: IPM parameter: tolerance

        ### References
        - "l1 Trend Filtering", S. Kim, K. Koh, ,S. Boyd and D. Gorinevsky https://web.stanford.edu/~gorin/papers/l1_trend_filter.pdf#page=13.55
        - https://web.stanford.edu/~boyd/l1_tf/
    """
    function l1tf(y, λ; α = 1e-2, β = 5e-1, μ = 2,
        MAXITER = 40, MAXLSITER = 20, atol = 1e-4, verbose = false)
        # Dimensions
        n = length(y); m = n - 2

        # Operator matrices
        D = spdiagm(m, n, 0 => ones(m), 1 => -2ones(m), 2 => ones(m))
        DDᵀ = D * D'
        Dy = D * y

        # Vectors
        z,  _z = zeros(m), zeros(m)
        μ1, _μ1 = ones(m), ones(m)
        μ2, _μ2 = ones(m), ones(m)
        f1, _f1 =  z .- λ,  z .- λ
        f2, _f2 = -z .- λ, -z .- λ

        # Variables
        t = 1e-10
        pobj = Inf
        dobj = 0
        step = Inf

        verbose && println("""--------------------------------------------
                    l1 trend filtering via primal-dual algorithm
                    Version 0.7 May 1 2007
                    Kwangmoo Koh, Seung-Jean Kim, Stephen Boyd
                    --------------------------------------------
                    Itr\tPrimalObj.\tDualObj.\tGap""")

        for itr in 0:MAXITER
            Dᵀz  = D'z
            DDᵀz = D * Dᵀz
            w    = Dy - (μ1 - μ2)

            # Why `dot` is not used: https://stackoverflow.com/a/70653944/12285249
            pobj1 = .5*w'*(DDᵀ\w) + λ*sum(μ1 + μ2)
            pobj2 = .5*Dᵀz'Dᵀz + λ*sum(abs.(Dy - DDᵀz))
            pobj = min(pobj1, pobj2)
            dobj = -.5*Dᵀz'*Dᵀz + Dy'*z
            gap = pobj - dobj

            verbose && @printf("%d\t%.4e\t%.4e\t%.4e\n", itr, pobj, dobj, gap)

            if gap ≤ atol
                return L1tf(y - D'z, true)
            end

            if step ≥ 0.2
                t = max(2 * m * μ / gap, 1.2 * t)
            end

            S = DDᵀ - spdiagm((μ1 ./ f1) + (μ2 ./ f2))
            r = -DDᵀz + Dy + (1/t) ./ f1 - (1/t) ./ f2
            dz = S \ r
            dμ1 = -(μ1 + ((1/t) .+ dz .* μ1) ./ f1)
            dμ2 = -(μ2 + ((1/t) .- dz .* μ2) ./ f2)

            resDual = DDᵀz - w
            resCent = [-μ1 .* f1 .- (1/t); -μ2 .* f2 .- (1/t)]
            residual = [resDual; resCent]

            negIdx1 = dμ1 .< 0
            negIdx2 = dμ2 .< 0
            step = 1
            if any(negIdx1)
                step = min(step, 0.99 * minimum(-μ1[negIdx1] ./ dμ1[negIdx1]))
            end
            if any(negIdx2)
                step = min(step, 0.99 * minimum(-μ2[negIdx2] ./ dμ2[negIdx2]))
            end

            for _ in 1:MAXLSITER
                _z  = z .+ step * dz
                _μ1 = μ1 .+ step .* dμ1
                _μ2 = μ2 .+ step .* dμ2
                _f1 =  _z .- λ
                _f2 = -_z .- λ

                _resDual = DDᵀ * _z - Dy + _μ1 - _μ2
                _resCent = [-_μ1 .* _f1 .- (1/t); -_μ2 .* _f2 .- (1/t)]
                _residual = [_resDual; _resCent]

                if max(maximum(_f1), maximum(_f2)) < 0 &&
                norm(_residual) ≤ (1 - α * step) * norm(residual)
                    break
                end
                step *= β
            end
            (z, μ1, μ2, f1, f2) = (_z, _μ1, _μ2, _f1, _f2)
        end

        # Fail to solve
        return L1tf(y - D'z, false)
    end

    include("snp500.jl")

    export L1tf, l1tf, snp500
end # module L1TrendFiltering

