# Keywords that describe how a file stores its image, not the image: dropped
# when it is copied into a plain FITS file, which writes its own.
const _STRUCTURAL_KEY = r"^(SIMPLE|XTENSION|BITPIX|NAXIS\d*|EXTEND|NEXTEND|PCOUNT|GCOUNT|TFIELDS|TTYPE\d+|TFORM\d+|TUNIT\d+|BSCALE|BZERO|BLANK|CHECKSUM|DATASUM|ZIMAGE|ZCMPTYPE|ZBITPIX|ZNAXIS\d*|ZTILE\d+|ZNAME\d+|ZVAL\d+|ZMASKCMP|ZSIMPLE|ZTENSION|ZEXTEND|ZBLOCKED|ZPCOUNT|ZGCOUNT|ZHECKSUM|ZDATASUM|ZQUANTIZ|ZDITHER0|ZBLANK)$"

# a header-only primary HDU can read as a 0 × 0 image
_is_image(hdu) = hdu isa ImageHDU && ndims(hdu) == 2 && all(>(0), size(hdu))

"""
    image_extensions(path) -> Vector{Int}

The HDUs of the FITS file at `path` that hold a 2-D image, by number (1 is
the primary): one for a plain or tile-compressed (`.fz`) frame, one per
amplifier or detector for a camera that writes several.
"""
image_extensions(path::AbstractString) = FITS(f -> [k for k in 1:length(f) if _is_image(f[k])], path, "r")

# An IRAF section, "[x1:x2,y1:y2]", as index ranges along FITS axes 1 and 2
# (reversed where x1 > x2); `nothing` if it is not one.
function _section(s)
    m = s isa AbstractString ? match(r"^\s*\[\s*(\d+):(\d+)\s*,\s*(\d+):(\d+)\s*\]\s*$", s) : nothing
    m === nothing && return nothing
    a, b, c, d = parse.(Int, m.captures)
    return (a:(a <= b ? 1 : -1):b, c:(c <= d ? 1 : -1):d)
end

_cards(h) = Tuple{String,Any,String}[(k, v, c) for (k, v, c) in zip(h.keys, h.values, h.comments) if !occursin(_STRUCTURAL_KEY, k)]

# `h`'s cards, then those of `primary` it lacks
function _merged_header(h, primary)
    own = _cards(h)
    names = Set(first.(own))
    cards = vcat(own, [c for c in _cards(primary) if !(c[1] in names) && !(c[1] in ("", "COMMENT", "HISTORY"))])
    return FITSHeader(first.(cards), getindex.(cards, 2), last.(cards))
end

function _physical(hdu)
    h = read_header(hdu)
    data = read(hdu)
    image = eltype(data) <: Union{Int8,UInt8,Int16,UInt16,Float32} ? Float32.(data) : Float64.(data)
    blank = haskey(h, "ZBLANK") ? h["ZBLANK"] : haskey(h, "BLANK") ? h["BLANK"] : nothing
    if blank !== nothing
        sentinel = blank * get(h, "BSCALE", 1.0) + get(h, "BZERO", 0.0)
        is_blank = image .== sentinel
        any(is_blank) && (image[is_blank] .= median(image[.!is_blank]))
    end
    return image, h
end

# One detector read out through several amplifiers, each extension's DATASEC
# placed at its DETSEC, less its overscan's median (BIASSEC) and, if every one
# has a GAIN, in electrons.
function _assembled(f, hdus)
    primary = read_header(f[1])
    parts = map(hdus) do k
        image, h = _physical(f[k])
        (image=image, header=h, data=_section(get(h, "DATASEC", nothing)), det=_section(get(h, "DETSEC", nothing)),
         bias=_section(get(h, "BIASSEC", nothing)), gain=Float64(get(h, "GAIN", 0.0)),
         binning=parse.(Int, split(string(get(h, "CCDSUM", get(primary, "CCDSUM", "1 1"))))))
    end
    any(p -> p.data === nothing || p.det === nothing, parts) && return nothing
    # DETSEC counts unbinned detector pixels, and a subarray's doesn't start at 1
    origin = [minimum(p -> minimum(p.det[a]), parts) for a in 1:2]
    binned(r, o, b) = fld(first(r) - o, b) + 1:step(r):fld(last(r) - o, b) + 1
    parts = [merge(p, (det=Tuple(binned(p.det[a], origin[a], p.binning[a]) for a in 1:2),)) for p in parts]
    n1 = maximum(p -> maximum(p.det[1]), parts); n2 = maximum(p -> maximum(p.det[2]), parts)
    image = fill(NaN32, n1, n2)
    electrons = all(p -> p.gain > 0, parts)
    span(r) = minimum(r):maximum(r)
    for p in parts
        level = p.bias === nothing ? 0.0 : median(p.image[span(p.bias[1]), span(p.bias[2])])
        image[p.det[1], p.det[2]] = (p.image[span(p.data[1]), span(p.data[2])] .- level) .* (electrons ? p.gain : 1.0)
    end
    header = _merged_header(parts[1].header, primary)
    drop = Set(["DATASEC", "DETSEC", "BIASSEC", "TRIMSEC", "CCDSEC", "AMPSEC", "EXTNAME", "EXTVER"])
    cards = Tuple{String,Any,String}[c for c in zip(header.keys, header.values, header.comments) if !(c[1] in drop) && c[1] != "GAIN"]
    electrons && push!(cards, ("GAIN", 1.0, "electrons: each amplifier's gain applied"))
    return image, FITSHeader(first.(cards), getindex.(cards, 2), last.(cards))
end

"""
    plain_frame(path; outdir, extension=nothing) -> String

The frame at `path` as a plain FITS file, its image in the primary HDU: what
every function here that takes frames' paths reads. A frame already plain is
returned as it is, or as [`clean_blank_pixels`](@ref) copies it if it has
`BLANK` pixels. Otherwise it is written to `outdir` as floating point,
`BLANK` pixels replaced by the median of the others:

- a tile-compressed (`.fz`) or multi-extension file's image: `extension`
  (an HDU number or `EXTNAME`), by default the first 2-D image;
- by default, a detector read out through several amplifiers, one image
  extension each with `DATASEC` and `DETSEC` (IRAF's convention), assembled:
  each amplifier's data section placed at its detector section (counted in
  unbinned pixels: `CCDSUM`'s binning and a subarray's origin are taken
  into account), less the median of its overscan (`BIASSEC`) and, if every amplifier has a `GAIN`,
  in electrons (`GAIN` then 1). What differs between amplifiers beyond
  that is the bias frames' and flat field's to remove
  ([`calibrate_frames`](@ref)).

The header is the image's own, without the keywords of its storage (a
compressed image's `ZNAXIS1`… become `NAXIS1`…), and with the primary
header's keywords it lacks, where cameras that write extensions keep the
time, the exposure and the pointing. An existing copy is reused.
[`analyze_set`](@ref) calls it on every frame.
"""
function plain_frame(path::AbstractString; outdir::AbstractString, extension=nothing)
    plain = extension === nothing && FITS(f -> _is_image(f[1]), path, "r")
    plain && return clean_blank_pixels(path; outdir)
    stem = replace(replace(basename(path), r"\.fz$"i => ""), r"\.(fits|fit|fts)$"i => "")
    out = joinpath(outdir, extension === nothing ? "$stem.fits" : "$stem-$extension.fits")
    isfile(out) && return out
    image, header = FITS(path, "r") do f
        hdus = [k for k in 1:length(f) if _is_image(f[k])]
        isempty(hdus) && throw(ArgumentError("$path: no 2-D image"))
        assembled = extension === nothing && length(hdus) > 1 ? _assembled(f, hdus) : nothing
        assembled === nothing || return assembled
        image, h = _physical(f[extension === nothing ? first(hdus) : extension])
        image, _merged_header(h, read_header(f[1]))
    end
    mkpath(outdir)
    FITS(f -> write(f, image; header), out, "w")
    return out
end
