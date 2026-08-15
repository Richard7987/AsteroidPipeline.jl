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
mkdir -p "$dest/science" "$dest/reference"

center="36.29821,2.05276"
base="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"

fetch() {
    local filefracday="$1" subdir="$2"
    local frac="${filefracday:8:6}" yr="${filefracday:0:4}" md="${filefracday:4:4}"
    local name="ztf_${filefracday}_000451_zr_c01_o_q1_sciimg.fits"
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
    fetch "$filefracday" science
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
    fetch "$filefracday" reference
done
