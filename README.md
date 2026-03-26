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
git clone <repository-url>
cd ships
pip install -r requirements.txt
```

## Usage

```bash
python source_highlight.py <catalog_file> <start_date> <end_date>
```

### Arguments

- `catalog_file`: Path to source catalog (FITS or TXT format). Supported file extensions: `.fits`, `.txt`
- `start_date`: Start date in YYYYMMDD format
- `end_date`: End date in YYYYMMDD format

### Example

```bash
python source_highlight.py ort_cat.fits 20260320 20260430
```

This generates an output file: `ships_ort_cat_20260320_20260430.txt`

## Output

The tool generates a formatted text file organized by observation date. For each date, it displays all visible sources with:

- **Source Name**: Identifier from the input catalog
- **RA (J2000 HMS)**: Right Ascension in Hours, Minutes, Seconds
- **Dec (J2000 DMS)**: Declination in Degrees, Arcminutes, Arcseconds
- **Elongation**: Angular separation from the Sun (degrees)
- **Rise Time**: Source rise time in IST
- **Set Time**: Source set time in IST

Sources within each date are sorted by minimum solar elongation (smallest first).

## Observation Site

The tool is configured for the Ooty Radio Telescope (ORT) in southern India:

- Longitude: 76.66° E
- Latitude: 11.38° N
- Height: 2240 m

## Input Catalog Requirements

Your source catalog must be in FITS (.fits) or TXT (.txt) format and include the following columns:

- `source_name`: Name/identifier of each source
- **RA Column** (one of): `ra`, `raj2000`, or `ra_j2000` (Right Ascension in J2000, HMS format)
- **Dec Column** (one of): `dec`, `decj2000`, or `dec_j2000` (Declination in J2000, DMS format)

**Note**: The tool searches for RA and Dec columns case-insensitively with multiple naming options for flexibility. Ensure your input catalog includes at least one of the acceptable column names for each coordinate.

## Technical Details

- Uses **Astropy** for coordinate transformations and time calculations
- **Coordinate Frames**: GCRS for solar elongation calculations, CIRS for rise/set times (accounts for precession and nutation)
- Transforms all source coordinates to GCRS and CIRS frames at observation dates for accurate calculations
- Converts between celestial (CIRS) and sidereal time systems
- Computes hour angles and rise/set times with proper atmospheric refraction considerations
- Implements solar elongation filtering for IPS-specific observation scheduling

## Author

Hardik Medhi