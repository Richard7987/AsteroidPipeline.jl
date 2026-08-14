#!/usr/bin/env bash
# Downloads the real ZTF frames used by real_data_demo.jl from IRSA's public
# archive (no authentication required). See real_data_demo.jl for details
# on why this particular field/night was chosen.
set -euo pipefail

dest="$(dirname "$0")/../data/real"
mkdir -p "$dest"

base="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"
for filefracday in 20210315117014 20210315119711 20210315122141; do
    frac="${filefracday:8:6}"
    name="ztf_${filefracday}_000451_zr_c01_o_q1_sciimg.fits"
    out="$dest/$name"
    if [[ -f "$out" ]]; then
        echo "already have $name"
        continue
    fi
    echo "fetching $name"
    curl -sf "$base/2021/0315/$frac/$name" -o "$out"
done
