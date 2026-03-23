import warnings
warnings.filterwarnings("ignore")

import numpy as np
import argparse
import os
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, get_sun, GCRS
from astropy.time import Time
import astropy.units as u

def get_units(data, col_name, default_unit):
    """Checks for existing unit metadata or infers it based on value ranges."""
    col = data[col_name]
    if col.unit is not None:
        return col.unit
    
    # Heuristic for RA: If any value > 24, it's almost certainly degrees.
    if col_name.lower() in ['ra', 'raj2000']:
        if np.nanmax(col) > 24:
            return u.deg
        return u.hourangle
    return default_unit

def run_ships(file_path, start_str, end_str):
    print(f"\n--- SHIPS: Source Highlighter for IPS ---\n")

    location = EarthLocation.from_geodetic(lon=76.66*u.deg, lat=11.38*u.deg, height=2240*u.m) #ORT Coordinates
    lat_rad = location.lat.rad
    base_name = os.path.splitext(os.path.basename(file_path))[0]

    try:
        start_time = Time(f"{start_str[:4]}-{start_str[4:6]}-{start_str[6:]}")
        end_time = Time(f"{end_str[:4]}-{end_str[4:6]}-{end_str[6:]}")
    except:
        print("Error: Dates must be YYYYMMDD.")
        return

    dates = start_time + np.arange((end_time - start_time).value + 1) * u.day

    if file_path.lower().endswith('.fits'):
        raw_data = Table.read(file_path, format='fits')
    else:
        raw_data = Table.read(file_path, format='ascii.guestimate')

    _, unique_indices = np.unique(raw_data['source_name'], return_index=True)
    data = raw_data[np.sort(unique_indices)]

    # Dynamic Unit Detection
    ra_unit = get_units(data, 'raj2000', u.hourangle)
    dec_unit = get_units(data, 'decj2000', u.deg)
    
    print(f"Detected units for {file_path}: RA = {ra_unit}, Dec = {dec_unit}")

    src_icrs = SkyCoord(ra=data['raj2000'], dec=data['decj2000'], 
                        unit=(ra_unit, dec_unit), frame='icrs')

    
    print(f"Processing {len(src_icrs)} sources for {len(dates)} days...")

    suns = get_sun(dates)
    seps = suns[:, None].separation(src_icrs[None, :]).deg
    mask = seps < 90
    date_indices, src_indices = np.where(mask)

    if len(date_indices) == 0:
        print("No matches found.")
        return

    obs_times = dates[date_indices]
    relevant_srcs = src_icrs[src_indices]
    relevant_suns = suns[date_indices]
    
    src_gcrs = relevant_srcs.transform_to(GCRS(obstime=obs_times))
    dra = (src_gcrs.ra.deg - relevant_suns.ra.deg + 360) % 360
    signs = np.where(dra < 180, 1, -1)
    signed_elongs = seps[date_indices, src_indices] * signs

    dec_rad = src_gcrs.dec.rad
    cos_H = (np.sin(0) - np.sin(lat_rad) * np.sin(dec_rad)) / (np.cos(lat_rad) * np.cos(dec_rad))
    cos_H = np.clip(cos_H, -1, 1)
    H_hours = np.degrees(np.arccos(cos_H)) / 15.0

    lst_mids = obs_times.sidereal_time('mean', longitude=location.lon).hour
    r_ist = (obs_times + ((src_gcrs.ra.hour - H_hours - lst_mids) % 24)*u.hour + 5.5*u.hour)
    s_ist = (obs_times + ((src_gcrs.ra.hour + H_hours - lst_mids) % 24)*u.hour + 5.5*u.hour)
    
    results_map = {}
    for i in range(len(date_indices)):
        idx = src_indices[i]
        s_name = data['source_name'][idx]
        
        # Format RA to HMS and Dec to DMS for the header
        ra_hms = src_icrs[idx].ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
        dec_dms = src_icrs[idx].dec.to_string(unit=u.deg, sep=':', precision=1, pad=True, alwayssign=True)
        
        if s_name not in results_map:
            results_map[s_name] = {
                'min_abs_elong': 999, 
                'entries': [], 
                'coords': f"RA (J2000 HMS): {ra_hms}  Dec (J2000 DMS): {dec_dms}"
            }
        
        elong = round(signed_elongs[i], 2)
        results_map[s_name]['entries'].append({
            'date': obs_times[i].to_value('iso', subfmt='date'),
            'elong': elong,
            'rise': r_ist[i].datetime.strftime('%H:%M'),
            'set': s_ist[i].datetime.strftime('%H:%M')
        })
        results_map[s_name]['min_abs_elong'] = min(results_map[s_name]['min_abs_elong'], abs(elong))

    sorted_sources = sorted(results_map.keys(), key=lambda x: results_map[x]['min_abs_elong'])

    out_name = f"ships_{base_name}_{start_str}_{end_str}.txt"
    with open(out_name, 'w') as f:
        f.write(f"SHIPS SCHEDULE | Catalog: {file_path}\n{'='*75}\n\n")
        for s_name in sorted_sources:
            info = results_map[s_name]
            f.write(f"SOURCE: {s_name:<20} | {info['coords']}\n")
            f.write(f"{'Date':<15} {'Elongation':<15} {'Rise (IST)':<15} {'Set (IST)':<15}\n")
            f.write("-" * 65 + "\n")
            for e in info['entries']:
                f.write(f"{e['date']:<15} {e['elong']:<15} {e['rise']:<15} {e['set']:<15}\n")
            f.write("\n" + "*"*75 + "\n\n")

    print(f"Success! Rise/Set times for {len(sorted_sources)} highlighted sources stored in {out_name}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SHIPS: Source Highlighter for IPS")
    parser.add_argument("file", help="Catalog path (FITS/Text)")
    parser.add_argument("start", help="Start YYYYMMDD")
    parser.add_argument("end", help="End YYYYMMDD")
    args = parser.parse_args()
    run_ships(args.file, args.start, args.end)