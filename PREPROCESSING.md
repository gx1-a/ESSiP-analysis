# Sentinel-1 and ERA5 preprocessing

## Sentinel-1

The analysis starts from Sentinel-1A/C Extra Wide (EW) swath Ground Range Detected (GRD) products.

The final SNAP graph used by the project is provided as `preprocessing/SNAP_preprocessing_graph.xml`. The only changes made to the archived graph are removal of machine-specific input/output paths and removal of the unused `externalDEMFile` path. The processing used SNAP's named **GETASSE30** DEM. The scientific processing parameters are otherwise unchanged.

The SNAP chain is:

1. Apply Orbit File — Sentinel Precise orbit, automatic download;
2. Thermal Noise Removal — HH/HV;
3. Calibration — linear Sigma0 for HH/HV (`outputImageScaleInDb = false`);
4. Refined Lee speckle filtering;
5. Terrain Correction — GETASSE30, bilinear DEM/image resampling, nominal spacing 39.8658876 m, local incidence angle saved;
6. geographic subset to approximately 160.0–162.8°E and 72.4–72.9°N;
7. GeoTIFF export.

The graph therefore contains the complete SNAP-to-GeoTIFF step; no additional undocumented scientific correction is required between SNAP and MATLAB.

For each source product, set the graph `Read` file to the downloaded Sentinel-1 GRD ZIP and the `Write` file to a GeoTIFF named

`Subset_<Sentinel-1-product-basename>_Orb_NR_Cal_Spk_TC.tif`.

Run `preprocessing/reproject_to_epsg3413_batch.m` on those GeoTIFFs. The script reads the GeoTIFF CRS, transforms the data to EPSG:3413, resamples to 40 m, and applies the scene-wise incidence-angle correction in the linear domain:

`sigma_corr = sigma * cos(theta) / cos(theta_ref)`

where `theta_ref` is the median valid local incidence angle for that scene. The corrected HH/HV fields are then converted to dB using `10*log10` and saved as processed MAT files.

The principal MAT variables are `HH_corr_dB`, `HV_corr_dB`, `xq`, `yq`, `xv`, and `yv`. Script 06 requires only `HH_corr_dB`, `xv`, and `yv` from the processed SAR scenes.

`01_derive_rough_ice_footprint.m` uses `statistics_sample_stride = 4` when estimating each scene's robust median/MAD normalization parameters. The persistence fields themselves are evaluated on the full valid image grid.

The exact SNAP software version used at the time of the original processing was not retained, but the operator chain and analysis-relevant settings are preserved in the supplied graph.

## ERA5

Scripts 05 and 08 use hourly ERA5 data covering 1 November 2025 through 30 April 2026 over the study region. The required instantaneous fields are `u10`, `v10`, `t2m`, `d2m`, `sp`, `msl`, and `siconc`. The project also uses the mean/accumulation-style fields `avg_sdlwrf`, `avg_sdswrf`, `avg_ishf`, `avg_slhtf`, `avg_snlwrf`, and `avg_snswrf`.

Place the ERA5 NetCDF files anywhere below `ERA5/`. Scripts 05 and 08 search recursively for files named `data_stream-oper_stepType-instant.nc` and `data_stream-oper_stepType-avg.nc`.
