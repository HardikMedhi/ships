# SHIPS: Source Highlighter for IPS

SHIPS is a Python tool for scheduling astronomical observations of celestial sources for Interplanetary Scintillation (IPS) studies. It reads source catalogs and computes rise/set times at the Ooty Radio Telescope (ORT), Tamil Nadu, India, with special emphasis on sources within 90° of the Sun—the optimal condition for IPS observations.

## Features

- **Flexible Catalog Input**: Supports both FITS and ASCII format source catalogs
- **IPS-Optimized Filtering**: Automatically identifies sources within 90° solar elongation, ideal for scintillation studies
- **IST Time Conversion**: Computes rise/set times in Indian Standard Time (IST) for ORT operations
- **Formatted Schedule**: Generates human-readable observation schedules with precise timing

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

- `catalog_file`: Path to source catalog (FITS or ASCII format). Must contain columns `raj2000`, `decj2000`, and `source_name`
- `start_date`: Start date in YYYYMMDD format
- `end_date`: End date in YYYYMMDD format

### Example

```bash
python source_highlight.py ort_cat.fits 20260320 20260430
```

This generates an output file: `ships_ort_cat_20260320_20260430.txt`

## Output

The tool generates a formatted text file containing:

- **Source Name**: Identifier from the input catalog
- **Coordinates**: RA (J2000 in HMS) and Dec (J2000 in DMS)
- **Date**: Observation date
- **Elongation**: Angular separation from the Sun (degrees)
- **Rise Time**: Source rise time in IST
- **Set Time**: Source set time in IST

Sources are sorted by minimum solar elongation, with sources having smaller minimum elongations listed first.

## Observation Site

The tool is configured for the Ooty Radio Telescope (ORT) in southern India:

- Longitude: 76.66° E
- Latitude: 11.38° N
- Height: 2240 m

## Input Catalog Requirements

Your source catalog must include the following columns:

- `source_name`: Name/identifier of each source
- `raj2000`: Right Ascension (J2000, in hours or degrees)
- `decj2000`: Declination (J2000, in degrees)

The tool will automatically detect coordinate units based on value ranges.

## Technical Details

- Uses **Astropy** for coordinate transformations and time calculations
- Converts between celestial (GCRS) and sidereal time systems
- Computes hour angles and rise/set times with proper atmospheric refraction considerations
- Implements solar elongation filtering for IPS-specific observation scheduling

## Author

Hardik Medhi (with the help of Gemini 3)