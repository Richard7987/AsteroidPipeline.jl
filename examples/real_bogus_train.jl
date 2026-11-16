#=
A real/false classifier for tracklets, trained on the tables that
real_bogus_dataset.jl writes and compared with the production filters:

    julia --project=. examples/real_bogus_train.jl

Declared before running:

Model. Logistic regression on the ten features, standardized with the
training sets' means and deviations, ridge penalty 1, fitted by Newton's
method. The flux ratio and the rate enter as log10 (the ratio capped at
1000); a missing axis ratio takes the training median. Injected and
catalogued asteroids are one class, real.

Validation. One whole set is left out at a time: nothing from its field
is seen in training. The threshold is the one that keeps, on the training
sets, as many reals as the production filters keep there; the count that
matters is then the false tracklets each lets through on the set left out,
with the reals each keeps beside it.
=#
using LinearAlgebra, Printf, Serialization, Statistics

const CACHE = joinpath(@__DIR__, "..", "data", "real", "real_bogus")
const NAMES = ["n", "snr_min", "snr_med", "log flux_ratio", "residual", "log rate", "masked", "peak_frac",
               "axis_ratio", "dipole"]

raw(r) = [r.n, r.snr_min, r.snr_med, log10(min(r.flux_ratio, 1000.0)), r.residual, log10(max(r.rate, 1e-3)),
          r.masked, r.peak_frac, r.axis_ratio, r.dipole]

function standardizer(X)
    med = [median(filter(isfinite, X[:, j])) for j in axes(X, 2)]
    F = copy(X)
    for j in axes(F, 2), i in axes(F, 1)
        isfinite(F[i, j]) || (F[i, j] = med[j])
    end
    μ = vec(mean(F; dims=1)); σ = max.(vec(std(F; dims=1)), 1e-9)
    X -> begin
        Z = copy(X)
        for j in axes(Z, 2), i in axes(Z, 1)
            isfinite(Z[i, j]) || (Z[i, j] = med[j])
        end
        hcat(ones(size(Z, 1)), (Z .- μ') ./ σ')
    end
end

sigmoid(z) = 1 / (1 + exp(-z))

function fit(Z, y; λ=1.0, iters=50)
    w = zeros(size(Z, 2))
    R = Diagonal(vcat(0.0, fill(λ, size(Z, 2) - 1)))
    for _ in 1:iters
        p = sigmoid.(Z * w)
        step = (Z' * (Z .* (p .* (1 .- p))) + R + 1e-9I) \ (Z' * (p - y) + R * w)
        w -= step
        norm(step) < 1e-8 && break
    end
    w
end

tables = [deserialize(joinpath(CACHE, f)) for f in sort(readdir(CACHE)) if endswith(f, ".jls")]
rows = reduce(vcat, [t.rows for t in tables if haskey(t, :rows)])
sets = unique(r.set for r in rows)
X = reduce(vcat, [raw(r)' for r in rows])
y = Float64[r.label != :false for r in rows]
prod = [r.production for r in rows]
setof = [r.set for r in rows]

@printf("%d tracklets in %d sets: %d real, %d false\n\n", length(rows), length(sets), count(==(1), y), count(==(0), y))
@printf("%-9s %5s %6s | %-13s | %-13s | %5s\n", "left out", "real", "false", "production", "classifier", "AUC")
@printf("%-9s %5s %6s | %5s %7s | %5s %7s |\n", "", "", "", "real", "false", "real", "false")
tot = zeros(Int, 6)
for s in sets
    te = setof .== s; tr = .!te
    st = standardizer(X[tr, :])
    w = fit(st(X[tr, :]), y[tr])
    ptr = sigmoid.(st(X[tr, :]) * w); pte = sigmoid.(st(X[te, :]) * w)
    keep = count(prod[tr] .& (y[tr] .== 1))
    thr = sort(ptr[y[tr].==1]; rev=true)[max(keep, 1)]
    real = y[te] .== 1; fals = .!real
    pos = pte[real]; neg = pte[fals]
    auc = isempty(pos) || isempty(neg) ? NaN : mean(p > q ? 1.0 : p == q ? 0.5 : 0.0 for p in pos, q in neg)
    c = [count(real), count(fals), count(prod[te] .& real), count(prod[te] .& fals),
         count((pte .>= thr) .& real), count((pte .>= thr) .& fals)]
    tot .+= c
    @printf("%-9s %5d %6d | %5d %7d | %5d %7d | %5.3f\n", s, c..., auc)
end
@printf("%-9s %5d %6d | %5d %7d | %5d %7d |\n", "all", tot...)

st = standardizer(X)
w = fit(st(X), y)
println("\nweights, all sets (standardized features):")
for (n, v) in zip(vcat("intercept", NAMES), w)
    @printf("  %-15s %+6.2f\n", n, v)
end
