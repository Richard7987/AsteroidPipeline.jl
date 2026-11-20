"""
    calibrate_frames(paths; outdir, bias=String[], dark=String[], flat=String[],
                     exptime_key="EXPTIME", extension=nothing) -> Vector{String}

Raw frames reduced with the calibration frames taken for them, written to
`outdir` (a plain FITS file each, its header with a `HISTORY` line) and
returned as paths, for [`analyze_set`](@ref). Any of `bias`, `dark` and
`flat` may be left out; each is combined by its pixels' median.

- master bias: the bias frames' median;
- master dark: the darks less the bias, per second of exposure
  (`exptime_key`), scaled to each frame's; without bias frames the darks are
  subtracted as they are, which suits darks of the frames' own exposure;
- master flat: the flats less bias and dark, each over its own median, then
  over their median's median. Pixels whose response is below a fifth of it,
  or not finite (dead pixels, a vignetted corner, an overscan strip), are set
  to `NaN`, which [`fill_value_mask`](@ref) masks.

Each frame becomes `(raw − bias − dark·t) / flat`. Overscan is not fitted
nor trimmed; frames, compressed or in extensions, are read as
[`plain_frame`](@ref) reads them (`extension` for a mosaic camera's).
"""
function calibrate_frames(paths::AbstractVector{<:AbstractString}; outdir::AbstractString,
                          bias::AbstractVector{<:AbstractString}=String[], dark::AbstractVector{<:AbstractString}=String[],
                          flat::AbstractVector{<:AbstractString}=String[], exptime_key::AbstractString="EXPTIME",
                          extension=nothing)
    scratch = mktempdir()
    try
        read_plain(p) = FITS(plain_frame(p; outdir=scratch, extension), "r") do f
            (image=Float64.(read(f[1])), header=read_header(f[1]))
        end
        shape = isempty(paths) ? nothing : size(read_plain(first(paths)).image)
        load(p) = let d = read_plain(p)
            shape === nothing || size(d.image) == shape ||
                throw(DimensionMismatch("$(basename(p)) is $(join(size(d.image), "x")) pixels and the frames $(join(shape, "x")): calibration frames must be read out as the frames are"))
            d
        end
        exptime(h) = haskey(h, exptime_key) ? Float64(h[exptime_key]) :
                     throw(ArgumentError("no $exptime_key in a frame's header: the dark can't be scaled"))
        combine(images) = median(cat(images...; dims=3); dims=3)[:, :, 1]
        master_bias = isempty(bias) ? nothing : combine([load(p).image for p in bias])
        scaled = master_bias !== nothing
        master_dark = isempty(dark) ? nothing : combine(map(dark) do p
            d = load(p)
            scaled ? (d.image .- master_bias) ./ exptime(d.header) : d.image
        end)
        remove(img, h) = let x = master_bias === nothing ? img : img .- master_bias
            master_dark === nothing ? x : x .- (scaled ? master_dark .* exptime(h) : master_dark)
        end
        master_flat = isempty(flat) ? nothing : let
            m = combine(map(flat) do p
                d = load(p)
                x = remove(d.image, d.header)
                x ./ median(filter(isfinite, x))
            end)
            m ./= median(filter(isfinite, m))
            m[.!isfinite.(m) .| (m .< 0.2)] .= NaN
            m
        end
        mkpath(outdir)
        return map(paths) do p
            d = load(p)
            image = remove(d.image, d.header)
            master_flat === nothing || (image ./= master_flat)
            used = join(filter(!isempty, [isempty(bias) ? "" : "$(length(bias)) bias",
                                          isempty(dark) ? "" : "$(length(dark)) dark",
                                          isempty(flat) ? "" : "$(length(flat)) flat"]), ", ")
            h = d.header
            cards = Tuple{String,Any,String}[(k, v, c) for (k, v, c) in zip(h.keys, h.values, h.comments) if !occursin(_STRUCTURAL_KEY, k)]
            push!(cards, ("HISTORY", nothing, "AsteroidPipeline calibrate_frames: $used"))
            header = FITSHeader(first.(cards), getindex.(cards, 2), last.(cards))
            out = joinpath(outdir, replace(replace(basename(p), r"\.fz$"i => ""), r"\.(fits|fit|fts)$"i => "") *
                                   (extension === nothing ? "" : "-$extension") * ".fits")
            FITS(f -> write(f, Float32.(image); header), out, "w")
            out
        end
    finally
        rm(scratch; recursive=true, force=true)
    end
end
