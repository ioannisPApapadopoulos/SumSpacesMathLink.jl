# This first computes an FFT approximation of the inverse Fourier Transform 
# and then adds corrections using Mathematica routines for the evaluation of
# the inverse Fourier Transform integral at the points xx1 for ̃U₋₁ and V₀ 
# and points xx2 for ̃U₀ and V₁ if stabilise is false and ̃Uₙ₊₁ and Vₙ₊₂ if stabilise is true. 
# These are then interpolated via Interpolations.jl

function fft_mathematica_supporter_functions(λ::Number, μ::Number, η::Number; W::Real=1000., δ::Real=0.001, 
    I::AbstractVector=[-1.,1.], N::Int=5, stabilise::Bool=false,
    xx1=-10:0.01:10, xx2=-1.05:0.01:1.05, maxrecursion=100)

    if λ == μ == η ≈ 0 error("Special case, use fft_supporter_functions instead.") end
    
    s = unique(2. ./ (I[2:end] - I[1:end-1]))
    if s != [1.0]
        error("mathematica_correction currently can only handle translations of the reference element [-1,1].")
    end
    
    (x, uS) = supporter_functions(λ, μ, η, W=W, δ=δ, s=s, N=N, stabilise=stabilise)
    (x1, x2, uS) = mathematica_corrections(λ, μ, η, x, uS, s, N, stabilise=stabilise, xx1=xx1, xx2=xx2, maxrecursion=maxrecursion)
    return interpolate_supporter_functions(x1, x2, uS, I, s)
end

# This function parses the output by Mathematica to something more useful for Julia.
function parse_mathematica(val)
    val1 = split(val.value,"`")
    val2 = split(val1[2],"^")
    val1 = parse(Float64, val1[1])
    if length(val2) > 1
        val1 = val1 * 10^parse(Float64,val2[2])
    end
    return val1
end

# This uses the Mathematica routine NIntegrate to compute the inverse Fourier Transform
# to high precision at the points xx1 for ̃U₋₁ and V₀ and points xx2 for ̃U₀ and V₁ if 
# stabilise is false and ̃Uₙ₊₁ and Vₙ₊₂ if stabilise is true. 

# As this is quite expensive, we use DelimitedFiles.jl to save the computed values as the
# function runs. 

# If you wish to run mathematica_correction without first running fft_supporter_functions
# then initiliaze as follows: uS = ([[0.]], [[0.]], [[0.]], [[0.]]); x = [0.0]; 

function mathematica_corrections(λ::Number, μ::Number, η::Number, x::AbstractVector, 
    uS::NTuple{4, Vector}, s::AbstractVector, N::Int; stabilise::Bool=false, 
    xx1::AbstractVector=-10:0.01:10, xx2::AbstractVector=-1.05:0.01:1.05, maxrecursion::Int=10)

    if s != [1.0]
        error("mathematica_correction currently can only handle translations of the reference element [-1,1].")
    end
    
    (yywT0, yyU_1, yywT1, yyU0) = uS
    xx1 = Array(xx1)
    xx2 = Array(xx2)
    x = Array(x)
    
    # Find uncommon points between x and xx1/xx2
    c1 = findall(!in(xx1), x)
    c2 = findall(!in(xx2), x)

    # Remove any common points from given support functions 
    ywT0 = [yywT0[1][c1]]; yU_1 = [yyU_1[1][c1]]; 
    ywT1 = [yywT1[1][c2]]; yU0 = [yyU0[1][c2]];
    x1 = x[c1]; x2 = x[c2]

    # Add extra points to be computed
    x1 = vcat(x1, xx1); perm1 = sortperm(x1); x1 = sort(x1); 
    x2 = vcat(x2, xx2); perm2 = sortperm(x2); x2 = sort(x2); 

    if ~isdir("uS-lmbda-$λ-mu-$μ-eta-$η")
        mkdir("uS-lmbda-$λ-mu-$μ-eta-$η")
    end

    if isfile("uS-lmbda-$λ-mu-$μ-eta-$η/uS-base.txt")
        tmp = readdlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-base.txt")
        ywT0[1] = tmp[1,:]; yU_1[1] = tmp[2,:];
    else
        for y in xx1
            print("Computing xx1 value: $y. \n")
            a1 = weval(W`Re[1/(2*Pi)*NIntegrate[Pi * BesselJ[0,k]/(λ + η*I*k - μ*I*Sign[k] + Abs[k]) * Exp[I y k], {k,-∞,∞}, WorkingPrecision -> 15, PrecisionGoal -> 12, MaxRecursion -> m]]`;λ=λ,η=η,μ=μ,y=y,m=maxrecursion)
            a1 = parse_mathematica(a1)

            ywT0 = [vcat(ywT0[1],[a1])]
    
            a2 = weval(W`Re[1/(2*Pi)*NIntegrate[I*Pi*k*BesselJ[0,Abs[k]]/Abs[k]/(λ + η*I*k - μ*I*Sign[k] + Abs[k]) * Exp[I y k], {k,-∞,∞}, WorkingPrecision -> 15, PrecisionGoal -> 12, MaxRecursion -> m]]`;λ=λ,η=η,μ=μ,y=y,m=maxrecursion)
            a2 = parse_mathematica(a2)
            yU_1 = [vcat(yU_1[1],[a2])]
        end
        ywT0[1] = ywT0[1][perm1]
        yU_1[1] = yU_1[1][perm1]
    end

    for y in xx2
        print("Computing xx2 value: $y. \n")
        if stabilise==true
            a3 = weval(W`Re[1/(2*Pi)*NIntegrate[(-I)^(N+2) * Pi * BesselJ[N+2,k]/(λ + η*I*k - μ*I*Sign[k] + Abs[k]) * Exp[I y k], {k,-∞,∞}, WorkingPrecision -> 15, PrecisionGoal -> 12, MaxRecursion -> m]]`;λ=λ,η=η,μ=μ,N=N,y=y,m=maxrecursion)
        else
            a3 = weval(W`Re[1/(2*Pi)*NIntegrate[- I *Pi * BesselJ[1,k]/(λ + η*I*k - μ*I*Sign[k] + Abs[k]) * Exp[I y k], {k,-∞,∞}, WorkingPrecision -> 15, PrecisionGoal -> 12, MaxRecursion -> m]]`;λ=λ,η=η,μ=μ,y=y,m=maxrecursion)
        end
        a3 = parse_mathematica(a3)
        ywT1 = [vcat(ywT1[1],[a3])]
        
        if stabilise==true
            a4 = weval(W`Re[1/(2*Pi)*NIntegrate[(-I)^(N+2) * Pi * BesselJ[N+2,k]/(-λ*I*Sign[k] - μ + η*Abs[k] - I*k) * Exp[I y k], {k,-∞,∞},  WorkingPrecision -> 15, PrecisionGoal -> 12, MaxRecursion -> m]]`;λ=λ,η=η,μ=μ,N=N,y=y,m=maxrecursion)
        else
            a4 = weval(W`Re[1/(2*Pi)*NIntegrate[- I *Pi * BesselJ[1,k]/(-λ*I*Sign[k] - μ + η*Abs[k] - I*k) * Exp[I y k], {k,-∞,∞}, WorkingPrecision -> 15, PrecisionGoal -> 12, MaxRecursion -> m]]`;λ=λ,η=η,μ=μ,y=y,m=maxrecursion)
        end
        a4 = parse_mathematica(a4)
        yU0 = [vcat(yU0[1],[a4])]

    end
    ywT1[1] = ywT1[1][perm2]
    yU0[1] = yU0[1][perm2]

    if ~isfile("uS-lmbda-$λ-mu-$μ-eta-$η/uS-base.txt")
        writedlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-base.txt", [real.(ywT0[1]), real.(yU_1[1])])
    end
    if stabilise==true
        writedlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-N-$N.txt", [x1, x2, real.(ywT0[1]), real.(yU_1[1]), real.(ywT1[1]), real.(yU0[1])])
    else
        writedlm("uS-lmbda-$λ-mu-$μ-eta-$η/uS-N-1.txt", [x1, x2, real.(ywT0[1]), real.(yU_1[1]), real.(ywT1[1]), real.(yU0[1])])
    end
    return (x1, x2, (ywT0, yU_1, ywT1, yU0))
end