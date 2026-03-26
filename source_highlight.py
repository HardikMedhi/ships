import warnings
warnings.filterwarnings("ignore")

import numpy as np
import argparse
import os
import csv
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, get_sun, GCRS, CIRS
from astropy.time import Time
import astropy.units as u

# Note: RA is assumed to be in Hours, Minutes, Seconds (HMS)
# and Dec is assumed to be in Degrees, Arcminutes, Arcseconds (DMS)

# ORT Location Constants
ORT_LOCATION = EarthLocation.from_geodetic(lon=76.66*u.deg, lat=11.38*u.deg, height=2240*u.m)

def parse_dates(start_str, end_str):
    """Parse date strings (YYYYMMDD format) into Time objects and generate date range."""
    try:
        start_time = Time(f"{start_str[:4]}-{start_str[4:6]}-{start_str[6:]}")
        end_time = Time(f"{end_str[:4]}-{end_str[4:6]}-{end_str[6:]}")
    except:
        print("Error: Dates must be YYYYMMDD.")
        return None
    
    dates = start_time + np.arange((end_time - start_time).value + 1) * u.day
    return dates

def load_catalog(file_path):
    """Load source catalog from FITS or TXT file."""
    file_ext = os.path.splitext(file_path)[1].lower()
    
    if file_ext == '.fits':
        try:
            return Table.read(file_path, format='fits')
        except Exception as e:
            print(f"Error: Unable to read FITS file. There might be a problem with the file type or extension.")
            print(f"Details: {e}")
            return None
    elif file_ext == '.txt':
        try:
            return Table.read(file_path, format='ascii.guestimate')
        except Exception as e:
            print(f"Error: Unable to read TXT file. There might be a problem with the file type or extension.")
            print(f"Details: {e}")
            return None
    else:
        print(f"Error: Unsupported file extension '{file_ext}'. Only .fits and .txt files are supported.")
        return None

def remove_duplicates(raw_data):
    """Remove duplicate sources from catalog, keeping first occurrence."""
    _, unique_indices = np.unique(raw_data['source_name'], return_index=True)
    return raw_data[np.sort(unique_indices)]

def find_coordinate_columns(data):
    """Find RA and Dec column names in catalog (case-insensitive)."""
    data_cols_lower = [col.lower() for col in data.colnames]
    
    ra_candidates = ['ra', 'raj2000', 'ra_j2000']
    dec_candidates = ['dec', 'decj2000', 'dec_j2000']
    
    ra_col = None
    dec_col = None
    
    for candidate in ra_candidates:
        if candidate in data_cols_lower:
            ra_col = data.colnames[data_cols_lower.index(candidate)]
            break
    
    for candidate in dec_candidates:
        if candidate in data_cols_lower:
            dec_col = data.colnames[data_cols_lower.index(candidate)]
            break
    
    if ra_col is None or dec_col is None:
        print(f"Error: Could not find RA/Dec columns in the catalog.")
        print(f"Available columns: {data.colnames}")
        print(f"Expected RA columns: {ra_candidates}")
        print(f"Expected Dec columns: {dec_candidates}")
        return None, None
    
    return ra_col, dec_col

def find_flux_column(data):
    """Find flux column in catalog (case-insensitive). Returns None if not found."""
    data_cols_lower = [col.lower() for col in data.colnames]
    flux_candidates = ['flux', 'flux_jy', 'flux_density', 's_peak']
    
    for candidate in flux_candidates:
        if candidate in data_cols_lower:
            flux_col = data.colnames[data_cols_lower.index(candidate)]
            print(f"Found flux column: {flux_col} (units: Jansky)")
            return flux_col
    
    print("Warning: Flux column not detected. Proceeding without flux information.")
    return None

def create_skycoord(data, ra_col, dec_col, dates):
    """Create SkyCoord objects in ICRS frame and transform to GCRS and CIRS at observation dates."""
    ra_unit = u.hourangle  # RA assumed to be in HMS
    dec_unit = u.deg       # Dec assumed to be in DMS
    
    src_icrs = SkyCoord(ra=data[ra_col], dec=data[dec_col], 
                        unit=(ra_unit, dec_unit), frame='icrs')
    
    # Transform to GCRS coordinates at all dates (for elongation calculations)
    src_gcrs_all = src_icrs.transform_to(GCRS(obstime=dates[:, None]))
    
    # Transform to CIRS coordinates at all dates (for rise/set time calculations with proper precession/nutation)
    src_cirs_all = src_icrs.transform_to(CIRS(obstime=dates[:, None]))
    
    return src_icrs, src_gcrs_all, src_cirs_all

def filter_by_solar_elongation(src_gcrs_all, dates, elongation_limit=90):
    """Filter sources by solar elongation criterion (within elongation_limit degrees of Sun)."""
    suns = get_sun(dates)
    seps = suns[:, None].separation(src_gcrs_all).deg
    mask = seps < elongation_limit
    date_indices, src_indices = np.where(mask)
    
    return date_indices, src_indices, seps

def calculate_rise_set_times(src_cirs, obs_times, location):
    """Calculate rise and set times in Indian Standard Time (IST) using CIRS frame."""
    lat_rad = location.lat.rad
    dec_rad = src_cirs.dec.rad
    
    cos_H = (np.sin(0) - np.sin(lat_rad) * np.sin(dec_rad)) / (np.cos(lat_rad) * np.cos(dec_rad))
    cos_H = np.clip(cos_H, -1, 1)
    H_hours = np.degrees(np.arccos(cos_H)) / 15.0
    
    lst_mids = obs_times.sidereal_time('mean', longitude=location.lon).hour
    r_ist = (obs_times + ((src_cirs.ra.hour - H_hours - lst_mids) % 24)*u.hour + 5.5*u.hour)
    s_ist = (obs_times + ((src_cirs.ra.hour + H_hours - lst_mids) % 24)*u.hour + 5.5*u.hour)
    
    return r_ist, s_ist

def calculate_elongations(src_gcrs, relevant_suns, seps):
    """Calculate signed solar elongations (positive: east, negative: west)."""
    dra = (src_gcrs.ra.deg - relevant_suns.ra.deg + 360) % 360
    signs = np.where(dra < 180, 1, -1)
    signed_elongs = seps * signs
    return signed_elongs

def build_results_by_date(date_indices, src_indices, data, src_icrs, obs_times, 
                          signed_elongs, r_ist, s_ist, flux_col=None):
    """Organize observation results by date."""
    results_by_date = {}
    source_coords_cache = {}
    
    for i in range(len(date_indices)):
        idx = src_indices[i]
        s_name = data['source_name'][idx]
        date_str = obs_times[i].to_value('iso', subfmt='date')
        
        # Cache source coordinates for first occurrence
        if s_name not in source_coords_cache:
            ra_hms = src_icrs[idx].ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
            dec_dms = src_icrs[idx].dec.to_string(unit=u.deg, sep=':', precision=1, pad=True, alwayssign=True)
            source_coords_cache[s_name] = f"RA (J2000 HMS): {ra_hms}  Dec (J2000 DMS): {dec_dms}"
        
        if date_str not in results_by_date:
            results_by_date[date_str] = []
        
        elong = round(signed_elongs[i], 2)
        entry = {
            'source_name': s_name,
            'elong': elong,
            'rise': r_ist[i].datetime.strftime('%H:%M'),
            'set': s_ist[i].datetime.strftime('%H:%M'),
            'coords': source_coords_cache[s_name]
        }
        
        # Add flux if available
        if flux_col is not None:
            entry['flux'] = float(data[flux_col][idx])
        
        results_by_date[date_str].append(entry)
    
    return results_by_date

def create_observation_schedule(results_by_date):
    """Create optimized observation schedule by binning sources by elongation and selecting highest flux per bin."""
    scheduled_by_date = {}
    
    for date_str, sources in results_by_date.items():
        scheduled_by_date[date_str] = []
        
        # Check if flux information is available
        if not sources or 'flux' not in sources[0]:
            return None  # Flux data required for scheduling
        
        # Bin sources by 5-degree elongation bins
        bins = {}
        for source in sources:
            # Find the bin (round to nearest 5 degrees)
            bin_key = round(source['elong'] / 5) * 5
            if bin_key not in bins:
                bins[bin_key] = []
            bins[bin_key].append(source)
        
        # Select highest flux source from each bin
        for bin_key in sorted(bins.keys()):
            bin_sources = bins[bin_key]
            # Select source with highest flux
            selected = max(bin_sources, key=lambda x: x['flux'])
            scheduled_by_date[date_str].append(selected)
    
    return scheduled_by_date

def write_output_file(out_name, file_path, results_by_date):
    """Write formatted observation schedule to output file."""
    sorted_dates = sorted(results_by_date.keys())
    
    with open(out_name, 'w') as f:
        f.write(f"SHIPS SCHEDULE | Catalog: {file_path}\n{'='*75}\n\n")
        for date_str in sorted_dates:
            f.write(f"DATE: {date_str}\n")
            f.write(f"{'Source Name':<25} {'RA (J2000 HMS)':<20} {'Dec (J2000 DMS)':<20} {'Elongation':<15} {'Rise (IST)':<15} {'Set (IST)':<15}\n")
            f.write("-" * 110 + "\n")
            
            # Sort sources by elongation within each day
            sorted_sources_for_date = sorted(results_by_date[date_str], key=lambda x: x['elong'])
            for entry in sorted_sources_for_date:
                coords_parts = entry['coords'].split('  ')
                ra_part = coords_parts[0].split(': ')[1] if len(coords_parts) > 0 else ''
                dec_part = coords_parts[1].split(': ')[1] if len(coords_parts) > 1 else ''
                f.write(f"{entry['source_name']:<25} {ra_part:<20} {dec_part:<20} {entry['elong']:<15} {entry['rise']:<15} {entry['set']:<15}\n")
            f.write("\n" + "*"*110 + "\n\n")
    
    return sorted_dates

def write_output_csv(out_name, results_by_date):
    """Write observation schedule to CSV file with columns: date, source_name, RA, Dec, rise, set, flux."""
    sorted_dates = sorted(results_by_date.keys())
    
    with open(out_name, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(['date', 'source_name', 'RA (J2000 HMS)', 'Dec (J2000 DMS)', 'rise (IST)', 'set (IST)', 'flux (Jy)'])
        
        # Write data sorted by date, then by elongation within each date
        for date_str in sorted_dates:
            sorted_sources_for_date = sorted(results_by_date[date_str], key=lambda x: x['elong'])
            for entry in sorted_sources_for_date:
                # Extract RA and Dec from coords string
                coords_parts = entry['coords'].split('  ')
                ra_part = coords_parts[0].split(': ')[1] if len(coords_parts) > 0 else ''
                dec_part = coords_parts[1].split(': ')[1] if len(coords_parts) > 1 else ''
                
                # Get flux if available, else NaN
                flux_val = entry.get('flux', np.nan)
                
                writer.writerow([
                    date_str,
                    entry['source_name'],
                    ra_part,
                    dec_part,
                    entry['rise'],
                    entry['set'],
                    flux_val
                ])
    
    return sorted_dates

def run_ships(file_path, start_str, end_str, schedule=False):
    """Main orchestration function for SHIPS scheduling."""
    print(f"\n--- SHIPS: Source Highlighter for IPS ---\n")

    # Parse dates
    dates = parse_dates(start_str, end_str)
    if dates is None:
        return
    
    # Load and validate catalog
    raw_data = load_catalog(file_path)
    if raw_data is None:
        return
    
    # Process catalog
    data = remove_duplicates(raw_data)
    ra_col, dec_col = find_coordinate_columns(data)
    if ra_col is None or dec_col is None:
        return
    
    print(f"Found RA column: {ra_col}, Dec column: {dec_col}")
    print(f"Using hardcoded units: RA = HMS, Dec = DMS")
    
    # Check for flux column
    flux_col = find_flux_column(data)
    
    # Create sky coordinates
    src_icrs, src_gcrs_all, src_cirs_all = create_skycoord(data, ra_col, dec_col, dates)
    print(f"Processing {len(src_icrs)} sources for {len(dates)} days...")
    
    # Filter by solar elongation
    date_indices, src_indices, seps = filter_by_solar_elongation(src_gcrs_all, dates)
    
    if len(date_indices) == 0:
        print("No matches found.")
        return
    
    # Extract relevant data for matched observations
    obs_times = dates[date_indices]
    src_gcrs = src_gcrs_all[date_indices, src_indices]
    src_cirs = src_cirs_all[date_indices, src_indices]
    relevant_suns = get_sun(dates)[date_indices]
    
    # Calculate elongations and rise/set times
    signed_elongs = calculate_elongations(src_gcrs, relevant_suns, seps[date_indices, src_indices])
    r_ist, s_ist = calculate_rise_set_times(src_cirs, obs_times, ORT_LOCATION)
    
    # Organize and output results
    results_by_date = build_results_by_date(date_indices, src_indices, data, src_icrs, 
                                             obs_times, signed_elongs, r_ist, s_ist, flux_col)
    
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Output as CSV file
    out_name_csv = f"all_sources_{base_name}_{start_str}_{end_str}.csv"
    write_output_csv(out_name_csv, results_by_date)
    print(f"Success! Schedule saved to {out_name_csv}\n")
    
    # Generate optimized observation schedule if requested
    if schedule:
        if flux_col is None:
            print("Warning: --schedule flag used but no flux column detected. Skipping schedule generation.")
        else:
            print("Generating optimized observation schedule...")
            scheduled = create_observation_schedule(results_by_date)
            if scheduled is None:
                print("Error: Could not create schedule without flux information.")
            else:
                schedule_out_name_csv = f"schedule_{base_name}_{start_str}_{end_str}.csv"
                write_output_csv(schedule_out_name_csv, scheduled)
                print(f"Optimized schedule (highest flux per 5° elongation bin) saved to {schedule_out_name_csv}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SHIPS: Source Highlighter for IPS")
    parser.add_argument("file", help="Catalog path (FITS/Text)")
    parser.add_argument("start", help="Start YYYYMMDD")
    parser.add_argument("end", help="End YYYYMMDD")
    parser.add_argument("--schedule", action="store_true", help="Generate optimized observation schedule (requires flux column)")
    args = parser.parse_args()
    run_ships(args.file, args.start, args.end, schedule=args.schedule)