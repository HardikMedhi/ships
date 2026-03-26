# SHIPS: Source Highlighter for IPS

SHIPS is a Python tool for scheduling astronomical observations of celestial sources for Interplanetary Scintillation (IPS) studies. It reads source catalogs and computes rise/set times at the Ooty Radio Telescope (ORT), Tamil Nadu, India, with special emphasis on sources within 90° of the Sun—the optimal condition for IPS observations.

## Features

- **Flexible Catalog Input**: Supports both FITS and TXT format source catalogs with flexible column naming
- **IPS-Optimized Filtering**: Automatically identifies sources within 90° solar elongation, ideal for scintillation studies
- **IST Time Conversion**: Computes rise/set times in Indian Standard Time (IST) for ORT operations
- **Date-Based Scheduling**: Generates schedules organized by observation date for convenient planning

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/HardikMedhi/ships.git
cd ships
pip install -r requirements.txt
```

## Usage

```bash
python source_highlight.py <catalog_file> <start_date> <end_date> [--schedule]
```

### Arguments

- `catalog_file`: Path to source catalog (FITS or TXT format). Supported file extensions: `.fits`, `.txt`
- `start_date`: Start date in YYYYMMDD format
- `end_date`: End date in YYYYMMDD format
- `--schedule` (optional): Generate an optimized observation schedule by selecting the highest-flux source in each 5° elongation bin. Requires a flux column in the catalog

### Examples

Generate complete observation schedule:
```bash
python source_highlight.py ort_cat.fits 20260320 20260430
```
This generates: `all_sources_ort_cat_20260320_20260430.csv`

Generate optimized schedule (highest flux per elongation bin):
```bash
python source_highlight.py ort_cat.fits 20260320 20260430 --schedule
```
This generates:
- `all_sources_ort_cat_20260320_20260430.csv` (all visible sources)
- `schedule_ort_cat_20260320_20260430.csv` (optimized schedule)

## Output

The tool generates CSV files organized by observation date with columns: `date`, `source_name`, `RA (J2000 HMS)`, `Dec (J2000 DMS)`, `rise (IST)`, `set (IST)`, `flux (Jy)`

Each row represents one source observation with:

- **date**: Observation date (YYYY-MM-DD)
- **source_name**: Identifier from the input catalog
- **RA (J2000 HMS)**: Right Ascension in Hours, Minutes, Seconds
- **Dec (J2000 DMS)**: Declination in Degrees, Arcminutes, Arcseconds
- **rise (IST)**: Source rise time in IST
- **set (IST)**: Source set time in IST
- **flux (Jy)**: Flux in Jansky (NaN if not available in catalog)

Data is organized by date (chronologically), then by elongation within each date (west to east).

### Flux Information

If a flux column is detected (looking for: `flux`, `flux_jy`, `flux_density`, or `s_peak`), it is loaded in memory with units assumed to be Jansky. This enables generation of optimized observation schedules.

### Optimized Schedule (--schedule flag)

When `--schedule` is used with flux data available, an additional file is generated that contains an optimized observation schedule. For each date:
- Sources are binned by elongation (5° bins)
- Only the highest-flux source is selected from each bin
- This minimizes observing conflicts while maximizing signal-to-noise

## Observation Site

The tool is configured for the Ooty Radio Telescope (ORT) in southern India:

- Longitude: 76.66° E
- Latitude: 11.38° N
- Height: 2240 m

## Input Catalog Requirements

Your source catalog must be in FITS (.fits) or TXT (.txt) format and include the following columns:

- `source_name`: Name/identifier of each source (required)
- **RA Column** (one of): `ra`, `raj2000`, or `ra_j2000` (Right Ascension in J2000, HMS format, required)
- **Dec Column** (one of): `dec`, `decj2000`, or `dec_j2000` (Declination in J2000, DMS format, required)
- **Flux Column** (one of, optional for scheduling): `flux`, `flux_jy`, `flux_density`, or `s_peak` (in Jansky units)

**Note**: The tool searches for RA and Dec columns case-insensitively with multiple naming options for flexibility. Flux information is automatically detected if present and is used for optimized scheduling with the `--schedule` flag.

## Technical Details

- Uses **Astropy** for coordinate transformations and time calculations
- **Coordinate Frames**: 
  - GCRS for solar elongation calculations (unaffected by Earth orientation)
  - CIRS for rise/set times (accounts for precession and nutation)
- Transforms all source coordinates to GCRS and CIRS frames at observation dates
- Converts between celestial (CIRS) and sidereal time systems
- Computes hour angles and rise/set times in IST (UTC+5:30)
- Implements solar elongation filtering for IPS-specific observation scheduling
- Optional flux-based scheduling: bins sources by 5° elongation intervals and selects highest flux per bin

## Author

Hardik Medhi