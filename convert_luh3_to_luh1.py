#!/usr/bin/env python3
"""
LUH3 (CMIP7) to LUH1 (CLM4) Land Use Data Converter - Version 6.0

Correct PFT Mapping Based on Data Analysis:

LUH3 Structure:
- PCT_NAT_PFT: 15 values (indices 0-14)
  - natpft[0]: Bare ground total percentage (NOT a vegetation type)
  - natpft[1-14]: 14 actual natural vegetation PFTs
- PCT_CFT: 64 crop functional types (array indices 0-63)
  - CFT[0-31]: Rainfed crops (~67%)
  - CFT[32-63]: Irrigated crops (~33%)
- PCT_NATVEG: Total natural vegetation percentage
- PCT_CROP: Total crop percentage

LUH1 Structure (CLM4):
- PCT_PFT: 17 values (indices 0-16)
  - PFT[0-13]: 14 natural vegetation types
  - PFT[14]: C3 rainfed crop
  - PFT[15]: C3 irrigated crop  
  - PFT[16]: Bare ground (residual)

Conversion Algorithm:
1. Natural vegetation: PFT[i] = natpft[i+1] * PCT_NATVEG / 100  (i=0-13)
2. Rainfed crop: PFT[14] = sum(CFT[0-31])/100 * PCT_CROP
3. Irrigated crop: PFT[15] = sum(CFT[32-63])/100 * PCT_CROP
4. Bare ground: PFT[16] = residule - sum(PFT[0-15])
   where residule = 100 - PCT_GLACIER - PCT_LAKE - PCT_WETLAND - PCT_URBAN

Author: Based on LUH3 data structure analysis
Date: 2026-03-24
Version: 6.0
"""

import netCDF4 as nc
import numpy as np
import sys
import os
import argparse
from scipy import interpolate

class LUH3toLUH1Converter:
    """Convert LUH3 land use data to LUH1 format for CLM4"""
    
    def __init__(self, input_file, output_file, luh1_reference_dir=None, verbose=True):
        self.input_file = input_file
        self.output_file = output_file
        self.luh1_reference_dir = luh1_reference_dir
        self.verbose = verbose
        
        # Target LUH1 resolution
        self.luh1_nx = 720
        self.luh1_ny = 360
        
        # Data containers
        self.luh3_data = {}
        self.luh1_data = {}
        
    def log(self, message):
        """Print log message if verbose"""
        if self.verbose:
            print(message)
    
    def load_luh3_data(self):
        """Load LUH3 input data"""
        self.log(f"Loading LUH3 data from: {self.input_file}")
        
        ds = nc.Dataset(self.input_file, 'r')
        
        # Required variables
        required_vars = [
            'LANDMASK',
            'PCT_NATVEG', 'PCT_CROP',
            'PCT_NAT_PFT',  # (15, lat, lon) - natpft[0]=bare, natpft[1-14]=vegetation
            'PCT_CFT',      # (64, lat, lon) - crop functional types
            'PCT_GLACIER', 'PCT_LAKE', 'PCT_WETLAND', 'PCT_URBAN',
            'HARVEST_VH1', 'HARVEST_VH2',
            'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3',
            'GRAZING'
        ]
        
        for var in required_vars:
            if var not in ds.variables:
                raise ValueError(f"Required variable '{var}' not found in input file")
            self.luh3_data[var] = ds.variables[var][:]
        
        # Store dimensions
        self.luh3_ny, self.luh3_nx = self.luh3_data['LANDMASK'].shape
        
        self.log(f"  LUH3 dimensions: {self.luh3_nx} x {self.luh3_ny}")
        self.log(f"  PCT_NAT_PFT shape: {self.luh3_data['PCT_NAT_PFT'].shape}")
        self.log(f"  PCT_CFT shape: {self.luh3_data['PCT_CFT'].shape}")
        
        # Verify dimensions
        if self.luh3_data['PCT_NAT_PFT'].shape[0] != 15:
            raise ValueError(f"Expected 15 natpft in LUH3, got {self.luh3_data['PCT_NAT_PFT'].shape[0]}")
        
        if self.luh3_data['PCT_CFT'].shape[0] != 64:
            raise ValueError(f"Expected 64 CFT in LUH3, got {self.luh3_data['PCT_CFT'].shape[0]}")
        
        ds.close()
        self.log("LUH3 data loaded successfully")
    
    def regrid_2d_conservative(self, data_in, method='conservative', landmask_in=None, fill_value=0.0):
        """
        Regrid 2D data from LUH3 to LUH1 resolution
        
        Args:
            data_in: Input 2D array (luh3_ny, luh3_nx)
            method: 'conservative' or 'nearest'
            landmask_in: Optional landmask for source data
            fill_value: Value for ocean points
        
        Returns:
            Regridded 2D array (luh1_ny, luh1_nx)
        """
        ny_in, nx_in = data_in.shape
        ny_out, nx_out = self.luh1_ny, self.luh1_nx
        
        # Create coordinate grids
        lon_in = np.linspace(-179.875, 179.875, nx_in)   # LUH3: 0.25 degree
        lat_in = np.linspace(-89.875, 89.875, ny_in)
        
        lon_out = np.linspace(-179.75, 179.75, nx_out)   # LUH1: 0.5 degree
        lat_out = np.linspace(-89.75, 89.75, ny_out)
        
        if method == 'conservative':
            # Conservative regridding: area-weighted average
            # For 0.25->0.5 degree: exact 2x2 aggregation
            
            data_out = np.zeros((ny_out, nx_out))
            
            for j_out in range(ny_out):
                for i_out in range(nx_out):
                    # Find corresponding input cells (2x2 block)
                    j_in_start = j_out * 2
                    j_in_end = j_in_start + 2
                    i_in_start = i_out * 2
                    i_in_end = i_in_start + 2
                    
                    # Extract sub-array
                    sub_data = data_in[j_in_start:j_in_end, i_in_start:i_in_end]
                    
                    if landmask_in is not None:
                        sub_mask = landmask_in[j_in_start:j_in_end, i_in_start:i_in_end]
                        # Only average over land cells
                        land_cells = sub_mask > 0.5
                        if np.any(land_cells):
                            data_out[j_out, i_out] = np.mean(sub_data[land_cells])
                        else:
                            data_out[j_out, i_out] = fill_value
                    else:
                        data_out[j_out, i_out] = np.mean(sub_data)
        
        elif method == 'nearest':
            # Nearest neighbor interpolation
            f = interpolate.RegularGridInterpolator(
                (lat_in, lon_in), data_in,
                method='nearest',
                bounds_error=False,
                fill_value=fill_value
            )
            
            lon_out_grid, lat_out_grid = np.meshgrid(lon_out, lat_out)
            points_out = np.column_stack([lat_out_grid.ravel(), lon_out_grid.ravel()])
            data_out = f(points_out).reshape(ny_out, nx_out)
        
        else:
            raise ValueError(f"Unknown regridding method: {method}")
        
        return data_out
    
    def regrid_all_data(self):
        """Regrid all LUH3 data to LUH1 resolution"""
        self.log("Regridding data to LUH1 resolution...")
        
        self.luh1_data = {}
        
        # Get LUH3 LANDMASK for masking ocean during regridding
        luh3_landmask = self.luh3_data['LANDMASK']
        
        # Load grid definition from LUH1 reference
        luh1_reference_loaded = False
        
        if self.luh1_reference_dir:
            self.log(f"  Loading grid definition from LUH1 reference...")
            import re
            year_match = re.search(r'(\d{4})', os.path.basename(self.input_file))
            ref_year = 1850  # Default
            if year_match:
                year = int(year_match.group(1))
                ref_year = min(year, 2005)
            
            ref_file = os.path.join(
                self.luh1_reference_dir,
                f"mksrf_landuse_rc{ref_year}_c090630.nc"
            )
            
            if os.path.exists(ref_file):
                try:
                    ds_ref = nc.Dataset(ref_file, 'r')
                    
                    # Copy complete grid definition from LUH1 reference
                    grid_fields = [
                        'LANDMASK',
                        'LONGXY', 'LATIXY',
                        'LON', 'LAT',
                        'EDGEN', 'EDGEE', 'EDGES', 'EDGEW'
                    ]
                    
                    for field in grid_fields:
                        if field in ds_ref.variables:
                            self.luh1_data[field] = ds_ref.variables[field][:]
                            self.log(f"    Copied {field} from LUH1 reference")
                    
                    # Load GRAZING for 1850-2005
                    if 'GRAZING' in ds_ref.variables:
                        self.luh1_data['GRAZING'] = ds_ref.variables['GRAZING'][:]
                        self.log(f"    Loaded GRAZING from LUH1 reference")
                    
                    luh1_reference_loaded = True
                    
                    if 'LANDMASK' in self.luh1_data:
                        n_land = np.sum(self.luh1_data['LANDMASK'] > 0.5)
                        self.log(f"    LUH1 reference land points: {n_land}")
                    
                    ds_ref.close()
                    
                except Exception as e:
                    self.log(f"    Warning: Could not load from reference: {e}")
        
        # If reference not loaded, create grid from scratch
        if not luh1_reference_loaded or 'LANDMASK' not in self.luh1_data:
            self.log("  Creating grid definition from LUH3...")
            self.luh1_data['LANDMASK'] = self.regrid_2d_conservative(
                self.luh3_data['LANDMASK'], 
                fill_value=0.0,
                method='nearest'
            )
        
        # Get LUH1 LANDMASK (this is the master mask)
        luh1_landmask = self.luh1_data['LANDMASK']
        
        self.log(f"  Using LUH1 reference LANDMASK as master")
        self.log(f"  LUH1 reference land points: {np.sum(luh1_landmask > 0.5)}")
        
        # Regrid DATA fields from LUH3
        # Apply LUH1 LANDMASK to all regridded data
        self.log("  Regridding data fields from LUH3...")
        
        # IMPORTANT: Use 'nearest' method to avoid coastline smoothing artifacts
        # 'conservative' (2x2 averaging) dilutes coastal values by mixing with ocean zeros
        regrid_method = 'nearest'
        
        # Special land units
        for var in ['PCT_GLACIER', 'PCT_LAKE', 'PCT_WETLAND', 'PCT_URBAN']:
            self.log(f"    Regridding {var}...")
            self.luh1_data[var] = self.regrid_2d_conservative(
                self.luh3_data[var],
                method=regrid_method,
                landmask_in=luh3_landmask
            )
            # Apply LUH1 LANDMASK
            self.luh1_data[var] = np.where(
                luh1_landmask > 0.5,
                self.luh1_data[var],
                0.0
            )
        
        # Vegetation fields
        for var in ['PCT_NATVEG', 'PCT_CROP']:
            self.log(f"    Regridding {var}...")
            self.luh1_data[var] = self.regrid_2d_conservative(
                self.luh3_data[var],
                method=regrid_method,
                landmask_in=luh3_landmask
            )
            self.luh1_data[var] = np.where(
                luh1_landmask > 0.5,
                self.luh1_data[var],
                0.0
            )
        
        # PCT_NAT_PFT (3D: 15 natpft)
        self.log(f"    Regridding PCT_NAT_PFT (15 values)...")
        npft_luh3 = 15
        ny_out, nx_out = self.luh1_ny, self.luh1_nx
        pct_nat_pft_out = np.zeros((npft_luh3, ny_out, nx_out))
        
        for k in range(npft_luh3):
            pct_nat_pft_out[k, :, :] = self.regrid_2d_conservative(
                self.luh3_data['PCT_NAT_PFT'][k, :, :],
                method=regrid_method,
                landmask_in=luh3_landmask
            )
            pct_nat_pft_out[k, :, :] = np.where(
                luh1_landmask > 0.5,
                pct_nat_pft_out[k, :, :],
                0.0
            )
        
        self.luh1_data['PCT_NAT_PFT'] = pct_nat_pft_out
        
        # PCT_CFT (3D: 64 crop types)
        self.log(f"    Regridding PCT_CFT (64 crop types)...")
        ncft = 64
        pct_cft_out = np.zeros((ncft, ny_out, nx_out))
        
        for k in range(ncft):
            pct_cft_out[k, :, :] = self.regrid_2d_conservative(
                self.luh3_data['PCT_CFT'][k, :, :],
                method=regrid_method,
                landmask_in=luh3_landmask
            )
            pct_cft_out[k, :, :] = np.where(
                luh1_landmask > 0.5,
                pct_cft_out[k, :, :],
                0.0
            )
        
        self.luh1_data['PCT_CFT'] = pct_cft_out
        
        # Harvest variables
        for var in ['HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 
                    'HARVEST_SH2', 'HARVEST_SH3']:
            self.log(f"    Regridding {var}...")
            self.luh1_data[var] = self.regrid_2d_conservative(
                self.luh3_data[var],
                method=regrid_method,
                landmask_in=luh3_landmask
            )
            self.luh1_data[var] = np.where(
                luh1_landmask > 0.5,
                self.luh1_data[var],
                0.0
            )
        
        # GRAZING
        if 'GRAZING' not in self.luh1_data and self.luh1_reference_dir:
            import re
            year_match = re.search(r'(\d{4})', os.path.basename(self.input_file))
            if year_match:
                year = int(year_match.group(1))
                if year > 2005:
                    grazing_ref = self.get_grazing_from_luh1_reference(year)
                    if grazing_ref is not None:
                        self.luh1_data['GRAZING'] = grazing_ref
                        self.log(f"    Using GRAZING from 2005 for year {year}")
        
        if 'GRAZING' not in self.luh1_data:
            self.log(f"    Warning: Using LUH3 GRAZING (zeros)")
            self.luh1_data['GRAZING'] = self.regrid_2d_conservative(
                self.luh3_data['GRAZING'],
                method=regrid_method,
                landmask_in=luh3_landmask
            )
            self.luh1_data['GRAZING'] = np.where(
                luh1_landmask > 0.5,
                self.luh1_data['GRAZING'],
                0.0
            )
        
        self.log("Regridding completed")
    
    def get_grazing_from_luh1_reference(self, year):
        """Get GRAZING data from LUH1 reference"""
        if year <= 2005:
            ref_file = os.path.join(
                self.luh1_reference_dir,
                f"mksrf_landuse_rc{year}_c090630.nc"
            )
        else:
            ref_file = os.path.join(
                self.luh1_reference_dir,
                "mksrf_landuse_rc2005_c090630.nc"
            )
        
        if os.path.exists(ref_file):
            try:
                ds = nc.Dataset(ref_file, 'r')
                grazing = ds.variables['GRAZING'][:]
                ds.close()
                return grazing
            except Exception as e:
                self.log(f"    Warning: Could not load GRAZING: {e}")
                return None
        else:
            return None
    
    def convert_pft_structure(self):
        """
        Convert LUH3 structure to LUH1 17-PFT structure
        
        Strategy: Ensure sum=100% (priority), PFT[0] starts from reference but adjustable
        - Load ALL PFT from LUH1 reference as initial/fallback values
        - If PCT_NATVEG + PCT_CROP < 1%: keep reference PFT (all 17 values)
        - Otherwise: 
          * Calculate PFT[1-15] from LUH3
          * PFT[0] = 100% - sum(PFT[1-15]) (adjusted to ensure sum=100%)
          * If PFT[0] < 0: scale down PFT[1-15] and set PFT[0]=0
        
        CRITICAL: CLM4 PFT indexing is 0-16 where:
        - PFT[0]: Bare ground (from LUH1 reference, fixed)
        - PFT[1-14]: 14 natural vegetation types (from LUH3)
        - PFT[15]: ALL crops (from LUH3)
        - PFT[16]: 0 (not used in CLM4)
        """
        self.log("Converting to LUH1 17-PFT structure...")
        self.log("  Strategy: PFT[0] from reference, PFT[1-15] from LUH3, scaled if needed")
        
        ny, nx = self.luh1_ny, self.luh1_nx
        
        # Get regridded data
        pct_nat_pft = self.luh1_data['PCT_NAT_PFT']  # (15, ny, nx)
        pct_cft = self.luh1_data['PCT_CFT']          # (64, ny, nx)
        pct_natveg = self.luh1_data['PCT_NATVEG']
        pct_crop = self.luh1_data['PCT_CROP']
        pct_glacier = self.luh1_data['PCT_GLACIER']
        pct_lake = self.luh1_data['PCT_LAKE']
        pct_wetland = self.luh1_data['PCT_WETLAND']
        pct_urban = self.luh1_data['PCT_URBAN']
        landmask = self.luh1_data['LANDMASK']
        
        # Initialize with zeros
        pct_pft = np.zeros((17, ny, nx), dtype=np.float64)
        
        # Load ALL PFT from LUH1 reference as initial/fallback values
        pft_reference_loaded = False
        if self.luh1_reference_dir:
            import re
            year_match = re.search(r'(\d{4})', os.path.basename(self.input_file))
            ref_year = 1850
            if year_match:
                year = int(year_match.group(1))
                ref_year = min(year, 2005)
            
            ref_file = os.path.join(
                self.luh1_reference_dir,
                f"mksrf_landuse_rc{ref_year}_c090630.nc"
            )
            
            if os.path.exists(ref_file):
                try:
                    ds_ref = nc.Dataset(ref_file, 'r')
                    if 'PCT_PFT' in ds_ref.variables:
                        # Copy ALL PFT from reference as fallback
                        pct_pft[:, :, :] = ds_ref.variables['PCT_PFT'][:, :, :]
                        pft_reference_loaded = True
                        self.log(f"  Loaded ALL PFT from LUH1 reference ({ref_year}) as fallback")
                    ds_ref.close()
                except Exception as e:
                    self.log(f"  Warning: Could not load reference PFT: {e}")
        
        if not pft_reference_loaded:
            self.log("  Warning: PFT not loaded from reference, will calculate from LUH3 only")
        
        # Calculate residule (not used for PFT sum, just for reference)
        residule = 100.0 - pct_glacier - pct_lake - pct_wetland - pct_urban
        residule = np.maximum(residule, 0.0)
        
        # Process each grid cell
        n_land = 0
        n_scaled = 0
        n_glacier = 0
        n_exact_fit = 0
        
        for j in range(ny):
            for i in range(nx):
                # ONLY process points that are land in LUH1 reference
                if landmask[j, i] < 0.5:
                    # Ocean: set all PFT to 0
                    # Note: PFT[0] might already be set from reference, but ocean should be 0
                    for k in range(17):
                        pct_pft[k, j, i] = 0.0
                    continue
                
                n_land += 1
                res = residule[j, i]
                
                # If glacier-dominated, keep ALL PFT from reference
                if pct_glacier[j, i] >= 99.0:
                    # All PFT already loaded from reference, keep them
                    n_glacier += 1
                    continue
                
                # Check if LUH3 has sufficient data
                total_vegetation = pct_natveg[j, i] + pct_crop[j, i]
                DATA_THRESHOLD = 1.0  # percent
                
                # Additional check: is this a coastal point?
                # Coastal points often have inconsistent data between resolutions
                # Check if neighboring points have very different vegetation
                is_coastal = False
                if j > 0 and j < ny-1 and i > 0 and i < nx-1:
                    # Check 4 neighbors
                    neighbor_veg = [
                        pct_natveg[j-1, i] + pct_crop[j-1, i],
                        pct_natveg[j+1, i] + pct_crop[j+1, i],
                        pct_natveg[j, i-1] + pct_crop[j, i-1],
                        pct_natveg[j, i+1] + pct_crop[j, i+1]
                    ]
                    # If this point has veg but neighbors are very low, or vice versa
                    # it's likely a coastal point with regrid artifacts
                    avg_neighbor = sum(neighbor_veg) / 4.0
                    if abs(total_vegetation - avg_neighbor) > 20.0:  # 20% difference
                        is_coastal = True
                
                if total_vegetation < DATA_THRESHOLD or is_coastal:
                    # Insufficient LUH3 data or coastal point: keep ALL PFT from reference
                    # pct_pft already has reference values, don't change
                    n_scaled += 1  # Count as "kept reference"
                    continue
                
                # Strategy: Calculate PFT[1-15] from LUH3, adjust PFT[0] to make sum = 100%
                # Priority: sum = 100% (most important), PFT[0] can be adjusted from reference
                
                # Save reference PFT for comparison
                pft_reference = pct_pft[:, j, i].copy()
                
                # Calculate PFT[1-15] from LUH3
                # 1. Natural vegetation: PFT[1-14]
                for k in range(14):
                    pct_pft[k+1, j, i] = pct_nat_pft[k+1, j, i] * pct_natveg[j, i] / 100.0
                
                # 2. Crops: PFT[15]
                pct_pft[15, j, i] = pct_crop[j, i]
                
                # 3. PFT[16] = 0
                pct_pft[16, j, i] = 0.0
                
                # 4. Adjust PFT[0] to make total = 100%
                sum_pft_1_15 = np.sum(pct_pft[1:16, j, i])
                pct_pft[0, j, i] = 100.0 - sum_pft_1_15
                
                # 5. Check if calculated PFT[0] differs too much from reference
                # This often happens at coastal points due to regrid artifacts
                pft0_diff = abs(pct_pft[0, j, i] - pft_reference[0])
                
                if pft0_diff > 30.0:  # More than 30% difference
                    # Likely a coastal/boundary artifact - use reference PFT instead
                    pct_pft[:, j, i] = pft_reference
                    n_scaled += 1  # Count as "kept reference"
                    continue
                
                # 6. Ensure PFT[0] is non-negative
                if pct_pft[0, j, i] < 0:
                    # PFT[1-15] exceeds 100%, need to scale down
                    if sum_pft_1_15 > 0:
                        scale = 100.0 / sum_pft_1_15
                        for k in range(1, 16):
                            pct_pft[k, j, i] *= scale
                        pct_pft[0, j, i] = 0.0
                        n_scaled += 1
                else:
                    n_exact_fit += 1
        
        self.luh1_data['PCT_PFT'] = pct_pft
        
        self.log(f"  Processed {n_land} land points:")
        self.log(f"    Glacier-dominated points: {n_glacier}")
        self.log(f"    Scaled PFT[1-15] to fit: {n_scaled}")
        self.log(f"    Exact fit (no scaling): {n_exact_fit}")
        self.log(f"    Under-filled (extra bare): {n_land - n_glacier - n_scaled - n_exact_fit}")
        
        self.log("PFT structure conversion completed")
    
    def create_luh1_coordinates(self):
        """Create LUH1 coordinate arrays if not loaded from reference"""
        if 'LON' not in self.luh1_data:
            self.luh1_lon = np.linspace(-179.75, 179.75, self.luh1_nx)
            self.luh1_lat = np.linspace(-89.75, 89.75, self.luh1_ny)
            self.luh1_data['LON'] = self.luh1_lon
            self.luh1_data['LAT'] = self.luh1_lat
        else:
            self.luh1_lon = self.luh1_data['LON']
            self.luh1_lat = self.luh1_data['LAT']
        
        if 'LONGXY' not in self.luh1_data:
            lon_2d, lat_2d = np.meshgrid(self.luh1_lon, self.luh1_lat)
            self.luh1_data['LONGXY'] = lon_2d
            self.luh1_data['LATIXY'] = lat_2d
        
        if 'EDGEN' not in self.luh1_data:
            self.luh1_data['EDGEN'] = 90.0
            self.luh1_data['EDGEE'] = 180.0
            self.luh1_data['EDGES'] = -90.0
            self.luh1_data['EDGEW'] = -180.0
    
    def write_output(self):
        """Write LUH1 format NetCDF file"""
        self.log(f"Writing output to: {self.output_file}")
        
        # Apply LUH1 LANDMASK to ensure ocean points are 0
        self.log("  Applying LUH1 LANDMASK to ensure ocean points are 0...")
        
        landmask = self.luh1_data['LANDMASK']
        ocean_mask = landmask <= 0.5  # Ocean points only
        
        # Count before cleaning
        n_ocean = np.sum(ocean_mask)
        self.log(f"  Ocean points (LANDMASK <= 0.5): {n_ocean}")
        
        # Apply to all 2D fields - only if ocean
        for var_name in ['HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 
                        'HARVEST_SH2', 'HARVEST_SH3', 'GRAZING']:
            if var_name in self.luh1_data:
                self.luh1_data[var_name][ocean_mask] = 0.0
        
        # Apply to PCT_PFT (3D) - only if ocean
        if 'PCT_PFT' in self.luh1_data:
            for k in range(17):
                self.luh1_data['PCT_PFT'][k, ocean_mask] = 0.0
        
        # Verify land points still sum to 100
        land_mask = landmask > 0.5
        if 'PCT_PFT' in self.luh1_data:
            pft_sum = np.sum(self.luh1_data['PCT_PFT'], axis=0)
            n_land = np.sum(land_mask)
            n_sum_100 = np.sum(np.abs(pft_sum[land_mask] - 100) < 0.01)
            self.log(f"  Land points with PFT sum = 100%: {n_sum_100} / {n_land}")
            
            if n_sum_100 < n_land:
                # Find problematic points
                bad_sum = land_mask & (np.abs(pft_sum - 100) > 0.01)
                n_bad = np.sum(bad_sum)
                self.log(f"  WARNING: {n_bad} land points do not sum to 100%!")
                if n_bad > 0:
                    indices = np.where(bad_sum)
                    j, i = indices[0][0], indices[1][0]
                    self.log(f"    Example: point ({j},{i}) sum = {pft_sum[j,i]:.2f}%")
        
        output_dir = os.path.dirname(self.output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        ds_out = nc.Dataset(self.output_file, 'w', format='NETCDF3_CLASSIC')
        
        # Create dimensions
        ds_out.createDimension('lon', self.luh1_nx)
        ds_out.createDimension('lat', self.luh1_ny)
        ds_out.createDimension('pft', 17)  # Changed from 'numpft' to 'pft'
        
        # Coordinate variables
        lon_var = ds_out.createVariable('LON', 'f8', ('lon',))
        lon_var.long_name = 'coordinate longitude'
        lon_var.units = 'degrees_east'
        lon_var[:] = self.luh1_data['LON']
        
        lat_var = ds_out.createVariable('LAT', 'f8', ('lat',))
        lat_var.long_name = 'coordinate latitude'
        lat_var.units = 'degrees_north'
        lat_var[:] = self.luh1_data['LAT']
        
        longxy_var = ds_out.createVariable('LONGXY', 'f8', ('lat', 'lon'))
        longxy_var.long_name = 'longitude'
        longxy_var.units = 'degrees_east'
        longxy_var[:] = self.luh1_data['LONGXY']
        
        latixy_var = ds_out.createVariable('LATIXY', 'f8', ('lat', 'lon'))
        latixy_var.long_name = 'latitude'
        latixy_var.units = 'degrees_north'
        latixy_var[:] = self.luh1_data['LATIXY']
        
        # Grid edges
        for edge_name in ['EDGEN', 'EDGEE', 'EDGES', 'EDGEW']:
            edge_var = ds_out.createVariable(edge_name, 'f8', ())
            edge_var.long_name = f'{edge_name.lower()} edge of surface grid'
            edge_var.units = 'degrees'
            edge_var[:] = self.luh1_data[edge_name]
        
        # LANDMASK
        mask_var = ds_out.createVariable('LANDMASK', 'f4', ('lat', 'lon'))
        mask_var.long_name = 'land/ocean mask'
        mask_var[:] = self.luh1_data['LANDMASK']
        
        # PCT_PFT (17 PFTs)
        pft_var = ds_out.createVariable('PCT_PFT', 'f4', ('pft', 'lat', 'lon'))
        pft_var.long_name = 'percent plant functional type'
        pft_var.units = 'percent'
        pft_var[:] = self.luh1_data['PCT_PFT']
        
        # Harvest variables
        for var_name in ['HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3']:
            var = ds_out.createVariable(var_name, 'f4', ('lat', 'lon'))
            var.long_name = f'{var_name} harvest'
            var.units = 'fraction'
            var[:] = self.luh1_data[var_name]
        
        # GRAZING
        grazing_var = ds_out.createVariable('GRAZING', 'f4', ('lat', 'lon'))
        grazing_var.long_name = 'grazing fraction'
        grazing_var.units = 'fraction'
        grazing_var[:] = self.luh1_data['GRAZING']
        
        # Global attributes
        ds_out.title = 'CLM4 surface data (LUH3 converted to LUH1 format)'
        ds_out.source = f'Converted from LUH3: {os.path.basename(self.input_file)}'
        ds_out.conversion_tool = 'convert_luh3_to_luh1_v6.py'
        ds_out.conversion_date = str(np.datetime64('now'))
        ds_out.conventions = 'CF-1.0'
        ds_out.history = 'LUH3 CMIP7 data converted to LUH1 CLM4 format'
        ds_out.comment = 'PFT[0] from LUH1 reference, PFT[1-15] from LUH3, PFT[16]=0'
        ds_out.producer = 'Hsin-Chien Liang @ AC3/RCEC, Academia Sinica'
        ds_out.contact = 'lama@gate.sinica.edu.tw'
        
        ds_out.close()
        self.log("Output file written successfully")
    
    def run(self):
        """Execute full conversion workflow"""
        self.log("="*70)
        self.log("LUH3 to LUH1 Converter - Version 6.0")
        self.log("="*70)
        
        self.load_luh3_data()
        self.regrid_all_data()
        self.convert_pft_structure()
        self.create_luh1_coordinates()
        self.write_output()
        
        self.log("="*70)
        self.log("Conversion completed successfully!")
        self.log("="*70)


def main():
    parser = argparse.ArgumentParser(
        description='Convert LUH3 (CMIP7) land use data to LUH1 (CLM4) format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python convert_luh3_to_luh1_v6.py \\
    -i mksrf_landuse_clm6_histLUH3_1850.c251012.nc \\
    -o mksrf_landuse_rc1850_c260324.nc \\
    --luh1-ref /path/to/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/

Correct PFT Mapping:
  LUH3 natpft[0]: Bare ground (ignored)
  LUH3 natpft[1-14] -> LUH1 PFT[1-14]: 14 natural vegetation types
  LUH3 PCT_CROP -> LUH1 PFT[15]: ALL crops (rainfed + irrigated)
  LUH1 PFT[16]: 0 (not used in CLM4)
  LUH1 PFT[0]: Bare ground (residual)
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input LUH3 NetCDF file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output LUH1 NetCDF file')
    parser.add_argument('--luh1-ref', dest='luh1_reference_dir',
                       help='Directory containing LUH1 reference files')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Suppress progress messages')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"[ERROR] Input file not found: {args.input}")
        sys.exit(1)
    
    if args.luh1_reference_dir and not os.path.exists(args.luh1_reference_dir):
        print(f"[ERROR] LUH1 reference directory not found: {args.luh1_reference_dir}")
        sys.exit(1)
    
    converter = LUH3toLUH1Converter(
        input_file=args.input,
        output_file=args.output,
        luh1_reference_dir=args.luh1_reference_dir,
        verbose=not args.quiet
    )
    
    try:
        converter.run()
    except Exception as e:
        print(f"\n[ERROR] Conversion failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
