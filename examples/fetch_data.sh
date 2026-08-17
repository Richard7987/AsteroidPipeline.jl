#!/usr/bin/env bash
# Downloads the real ZTF frames used by real_data_demo.jl from IRSA's public
# archive (no authentication required). See real_data_demo.jl for details on
# why this field/night/reference set was chosen.
#
# All frames are fetched as server-side cutouts (?center=...&size=1024pix)
# centred on the same sky position, so science and reference frames land on
# one shared ~1025x1025 footprint without any extra alignment step here
# (per-frame dithering is still ~5-6 px and is handled by reprojection in
# src/reference.jl, not by this script).
set -euo pipefail

dest="$(dirname "$0")/../data/real"
mkdir -p "$dest/science" "$dest/reference" "$dest/variable"

base="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"

fetch() {
    local filefracday="$1" subdir="$2" center="$3" field="$4" ccdid="$5" qid="$6" filtercode="$7"
    local frac="${filefracday:8:6}" yr="${filefracday:0:4}" md="${filefracday:4:4}"
    local name="ztf_${filefracday}_$(printf '%06d' "$field")_${filtercode}_c$(printf '%02d' "$ccdid")_o_q${qid}_sciimg.fits"
    local out="$dest/$subdir/$name"
    # Full-size cutouts are ~4 MB; a much smaller file here means a prior
    # run was interrupted mid-download (this happened once with a plain
    # connection reset) and left a truncated FITS file that looks present
    # but fails to open.
    if [[ -f "$out" ]] && [[ $(stat -c%s "$out" 2>/dev/null || stat -f%z "$out") -gt 1000000 ]]; then
        echo "already have $subdir/$name"
        return
    fi
    echo "fetching $subdir/$name"
    curl -sf "$base/$yr/$md/$frac/$name?center=$center&size=1024pix&gzip=false" -o "$out"
}

# Science: 2019-10-23, field 451/CCD 1/quadrant 1/zr, 5 exposures spread
# across the night's 6.25 h span (seeing 1.8-2.0 px, maglim 20.1-20.5).
# Contains known asteroid 135992 "2002 UY45" (Mv 18.4, ~1.9 mag above this
# night's detection limit) moving ~222" across the frames.
for filefracday in \
    20191023195590 20191023267002 20191023347164 20191023425208 20191023455822
do
    fetch "$filefracday" science "36.29821,2.05276" 451 1 1 zr
done

# Reference: 30 frames from other nights (2018-09 onward — earlier ZTF
# frames are listed in the archive's metadata but return 404 on download),
# best seeing first, for a median-stacked deep reference.
for filefracday in \
    20190805440637 20251222148287 20191116296343 20190904419977 20201018360370 \
    20200920406505 20210209152593 20200825428738 20241025416968 20190804500023 \
    20191009420243 20210810465208 20200128118102 20210904465081 20200816486516 \
    20210904466019 20220810478565 20210108206829 20251112309039 20241014277454 \
    20201023337419 20190804485162 20210913418067 20200915420046 20191116312153 \
    20251031279028 20201022416482 20250912391227 20250924338356 20210902478102
do
    fetch "$filefracday" reference "36.29821,2.05276" 451 1 1 zr
done

# Variable-star validation (see examples/variable_star_demo.jl): field
# 487/CCD 12/quadrant 1/zr, night 2019-06-10, a real high-cadence ZTF
# campaign with 144 exposures spanning 2.45 h. Thinned to every 5th
# exposure (~29 frames) — this raw-path demo has no reference stack to
# build, so many more frames aren't the bottleneck they are for
# real_data_demo.jl, but 144 near-duplicate 30 s exposures of the same
# ~2.45 h window add little beyond what 29 spread across it already give.
# Centred on ASASSN-V J183620.31 (VSX: type EW, period 0.322 d, amplitude
# 0.63 mag), the real variable this demo validates against.
for filefracday in \
    20190610318264 20190610358160 20190610360417 20190610362674 20190610364942 \
    20190610367199 20190610369456 20190610371713 20190610373958 20190610376215 \
    20190610378472 20190610380729 20190610382998 20190610385255 20190610387512 \
    20190610389768 20190610392025 20190610394282 20190610396539 20190610398796 \
    20190610401053 20190610403310 20190610405567 20190610407824 20190610410081 \
    20190610412326 20190610414606 20190610416863 20190610419120
do
    fetch "$filefracday" variable "279.085,5.56101" 487 12 1 zr
done
