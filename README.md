# LUH3 to LUH1 Converter for TaiESM1/CLM4

Convert LUH3 (CMIP7) land use data to LUH1 format for use with the TaiESM1/CLM4 `mksurfdat` tool.

**Version:** 6.0 | **Author:** Hsin-Chien Liang @ AC3/RCEC, Academia Sinica | **Contact:** lama@gate.sinica.edu.tw

---

## Overview

This tool reads LUH3-format NetCDF files (0.25°, CMIP7) and produces LUH1-format NetCDF files (0.5°, CLM4-compatible) suitable for driving TaiESM1/CESM1 land surface initialization and dynamic land use runs via `mksurfdat`.

---

## Data Format Comparison

| Property | LUH3 (Input) | LUH1 (Output) |
|---|---|---|
| Resolution | 1440×720 (0.25°) | 720×360 (0.5°) |
| PFT structure | 15 natural PFTs + 64 CFTs (separate) | 17 PFTs (combined) |
| Data type | double | float |
| FillValue | −9999.0 | −999.0 |

---

## PFT Mapping

### LUH3 Structure
- `PCT_NAT_PFT[0]`: bare ground total percentage (not a vegetation type)
- `PCT_NAT_PFT[1–14]`: 14 natural vegetation PFTs
- `PCT_CFT[0–31]`: rainfed crop functional types (~67%)
- `PCT_CFT[32–63]`: irrigated crop functional types (~33%)

### LUH1 Output Structure
| LUH1 Index | Source | Description |
|---|---|---|
| PFT[0] | LUH1 reference (adjusted) | Bare ground (residual) |
| PFT[1–14] | `PCT_NAT_PFT[1–14] × PCT_NATVEG / 100` | 14 natural vegetation types |
| PFT[15] | `PCT_CROP` | All crops (rainfed + irrigated combined) |
| PFT[16] | 0 (fixed) | Unused in CLM4 |

### Conversion Algorithm
```
1. residule = 100 - PCT_GLACIER - PCT_LAKE - PCT_WETLAND - PCT_URBAN

2. PFT[1–14] (natural vegetation):
   PCT_PFT[k] = PCT_NAT_PFT[k+1] × PCT_NATVEG / 100   (k = 0..13)

3. PFT[15] (crops):
   PCT_PFT[15] = PCT_CROP

4. PFT[16] = 0

5. PFT[0] (bare ground):
   PCT_PFT[0] = 100 - sum(PFT[1:16])

6. If PFT[0] < 0 (sum exceeds 100%):
   scale PFT[1–15] by 100 / sum(PFT[1:16]), set PFT[0] = 0

7. Coastal / low-vegetation points (PCT_NATVEG + PCT_CROP < 1%):
   keep all PFT values from LUH1 reference
```

---

## Spatial Regridding

Uses **nearest-neighbour** regridding (not area-averaging) to avoid coastal dilution artifacts when going from 0.25° to 0.5°. Each output grid point maps to the single nearest input point.

If a LUH1 reference directory is provided, the LANDMASK is taken directly from the reference files and used as the master mask throughout.

---

## Installation

```bash
# Python 3.6+
pip install numpy netCDF4 scipy
```

---

## Usage

### Single file
```bash
python convert_luh3_to_luh1.py \
    -i mksrf_landuse_clm6_histLUH3_1850.c251012.nc \
    -o mksrf_landuse_rc1850_c260324.nc \
    --luh1-ref /path/to/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/
```

### Full argument reference
```
-i / --input       Input LUH3 NetCDF file (required)
-o / --output      Output LUH1 NetCDF file (required)
--luh1-ref         Directory containing LUH1 reference files
                   (e.g. mksrf_landuse_rc1850_c090630.nc)
                   Used for: LANDMASK, grid definition, PFT[0] baseline, GRAZING
-q / --quiet       Suppress progress messages
```

### Batch conversion (1850–2023)

Edit `batch_convert_luh3_to_luh1.sh` to set input/output paths and file naming pattern, then run:

```bash
chmod +x batch_convert_luh3_to_luh1.sh
./batch_convert_luh3_to_luh1.sh
```

---

## LUH1 Reference Files

The `--luh1-ref` directory is strongly recommended. The converter uses it for:

- **LANDMASK** — master land/ocean mask applied to all output fields
- **Grid definition** — `LONGXY`, `LATIXY`, `EDGEN/E/S/W`
- **PFT[0] baseline** — bare ground fraction for glacier-dominated and data-sparse points
- **GRAZING** — for years ≤ 2005; for years > 2005, the 2005 reference value is reused

Expected filename pattern: `mksrf_landuse_rc{YEAR}_c090630.nc`

---

## Output Variables

| Variable | Dimensions | Description |
|---|---|---|
| `PCT_PFT` | (17, lat, lon) | Percent plant functional type |
| `HARVEST_VH1/VH2` | (lat, lon) | Harvest: virgin/heavily wooded |
| `HARVEST_SH1/SH2/SH3` | (lat, lon) | Harvest: secondary wooded/non-wooded |
| `GRAZING` | (lat, lon) | Grazing fraction |
| `LANDMASK` | (lat, lon) | Land/ocean mask |
| `LONGXY`, `LATIXY` | (lat, lon) | 2D coordinate arrays |
| `LON`, `LAT` | (lon), (lat) | 1D coordinate arrays |

---

## Output Verification

```bash
# 1. Check dimensions
ncdump -h output.nc | grep dimensions
# Expected: lon = 720, lat = 360, pft = 17

# 2. Check PFT sum (land points should be ~100%)
ncap2 -s 'pft_sum=PCT_PFT.total($pft)' output.nc test.nc

# 3. Check data range
ncdump -v PCT_PFT output.nc | grep -A 5 "PCT_PFT ="
# All values should be in 0–100
```

---

## Known Issues

### 1. Harvest units
- LUH3 uses `gC/m²/yr`; LUH1 expects unitless fractions
- Current version retains original values without unit conversion
- May require a scaling factor depending on what `mksurfdat` expects

### 2. 2005–2006 transition
- If original LUH1 data is available for 1850–2005, it is recommended to:
  - 1850–2005: use original LUH1 files
  - 2006–2023: use converted LUH3 files
  - Check continuity at the boundary year

### 3. Coastal/boundary artifacts
- Points where `|PCT_PFT[0]_calculated - PCT_PFT[0]_reference| > 30%` are flagged as coastal artifacts and reverted to reference PFT values
- High-latitude (polar) regions with large glacier/lake fractions are handled separately

### 4. Memory
- Single year: ~500 MB RAM
- Batch run: 2 GB+ recommended

---

## Next Steps: mksurfdat

```bash
mksurfdat \
    --luh1-dir /path/to/converted/luh1/ \
    --start-year 1850 \
    --end-year 2023 \
    --output fsurdat.nc \
    --output-dyndat fdyndat.nc
```

Refer to TaiESM1/CLM4 documentation for the full `mksurfdat` configuration.

---

## Troubleshooting

| Error | Likely Cause | Fix |
|---|---|---|
| `Required variable 'X' not found` | Missing variable in LUH3 input | Check input file with `ncdump -h` |
| `PCT_NAT_PFT shape[0] != 15` | Wrong LUH3 format version | Verify data source |
| `PCT_PFT sum > 100%` | Normal intermediate state | Script auto-rescales; check `residule` if persists in output |
| Out of memory | Large batch job | Reduce parallelism or add swap |

---

## References

1. R script — original CLM4 PFT mapping logic
2. CTSM5.2 land use data tool — percentage scaling logic (C)
3. LUH3 documentation: https://luh.umd.edu/
4. TaiESM1: Taiwan Earth System Model (AC3/RCEC, Academia Sinica)

---

## Version History

| Version | Date | Changes |
|---|---|---|
| 6.0 | 2026-03-24 | Corrected PFT mapping based on LUH3 data structure analysis; switched to nearest-neighbour regridding; LUH1 reference LANDMASK as master; coastal artifact detection |
| 1.0 | 2026-01-28 | Initial release — basic LUH3→LUH1 conversion, PFT structure conversion, spatial regridding, harvest data |

---

## License

Developed under the TaiESM1 project for scientific research purposes.

