%% reproject_to_epsg3413_batch.m
% Purpose: Reproject three-band GeoTIFFs to EPSG:3413,
% resample them to a 40 m grid, apply scene-wise cosine incidence angle
% normalization in the linear domain, and save the processed MAT files.
% Input bands: 1 = linear HH, 2 = linear HV, 3 = local incidence angle.
% The source GeoTIFF may be projected (as written by the supplied SNAP graph)
% or geographic; its CRS is read from the GeoTIFF metadata.
% MATLAB: R2022b

project_root = fileparts(fileparts(mfilename('fullpath')));
input_dir = fullfile(project_root,'input_geotiff');
output_dir = fullfile(project_root,'processed');

if ~isfolder(input_dir)
    error('Input GeoTIFF directory not found: %s',input_dir);
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

target_res = 40;       % target grid spacing (m)
proj3413   = projcrs(3413);

file_list = dir(fullfile(input_dir, '*.tif'));
nFiles    = numel(file_list);
fprintf('Found %d input scenes.\n', nFiles);

fail_log = {};   % failed input files

for k = 1 : nFiles
    fname = file_list(k).name;
    fpath = fullfile(input_dir, fname);
    fprintf('[%d/%d] %s\n', k, nFiles, fname);

    %% Read the three-band GeoTIFF
    try
        raw  = readgeoraster(fpath);
        info = georasterinfo(fpath);
    catch ME
        fprintf('    Read failed; skipping scene: %s\n', ME.message);
        fail_log{end+1} = fname;
        continue
    end

    R_src   = info.RasterReference;
    HH_raw  = double(raw(:,:,1));
    HV_raw  = double(raw(:,:,2));
    ANG_raw = double(raw(:,:,3));
    clear raw

    HH_raw(HH_raw   <= 0) = NaN;
    HV_raw(HV_raw   <= 0) = NaN;
    ANG_raw(ANG_raw <= 0) = NaN;

    %% Scene-median reference incidence angle
    theta_ref = nanmedian(ANG_raw(:));
    fprintf('    theta_ref = %.2f deg\n', theta_ref);

    %% Construct source pixel-center coordinates and transform to EPSG:3413
    [nR, nC] = size(HH_raw);
    col_idx  = 1:nC;
    row_idx  = (1:nR)';

    if isa(R_src, 'map.rasterref.MapCellsReference')
        x_native = R_src.XWorldLimits(1) + ...
            (col_idx - 0.5) * R_src.CellExtentInWorldX;
        y_native = R_src.YWorldLimits(2) - ...
            (row_idx - 0.5) * R_src.CellExtentInWorldY;
        [X_native, Y_native] = meshgrid(x_native,y_native);

        src_crs = info.CoordinateReferenceSystem;
        if isempty(src_crs) || ~isa(src_crs,'projcrs')
            error(['Projected input lacks a usable projected CRS in its ' ...
                'GeoTIFF metadata: %s'],fname);
        end
        [LAT_s, LON_s] = projinv(src_crs,X_native,Y_native);
        clear X_native Y_native

    elseif isa(R_src, 'map.rasterref.GeographicCellsReference')
        lon_native = R_src.LongitudeLimits(1) + ...
            (col_idx - 0.5) * R_src.CellExtentInLongitude;
        lat_native = R_src.LatitudeLimits(2) - ...
            (row_idx - 0.5) * R_src.CellExtentInLatitude;
        [LON_s, LAT_s] = meshgrid(lon_native,lat_native);

    else
        fprintf('    Unsupported raster-reference type; skipping scene.\n');
        fail_log{end+1} = fname;
        continue
    end

    [xs, ys] = projfwd(proj3413,LAT_s,LON_s);
    clear LAT_s LON_s

    %% Create a 40 m target grid in EPSG:3413
    valid_xy = isfinite(xs) & isfinite(ys);
    if ~any(valid_xy(:))
        fprintf('    No valid projected source coordinates; skipping scene.\n');
        fail_log{end+1} = fname;
        continue
    end

    xv = min(xs(valid_xy)) : target_res : max(xs(valid_xy));
    yv = (max(ys(valid_xy)) : -target_res : min(ys(valid_xy)))';
    [xq, yq] = meshgrid(xv, yv);

    src_valid = valid_xy & isfinite(HH_raw) & isfinite(HV_raw) & isfinite(ANG_raw);
    if nnz(src_valid) < 4
        fprintf('    Too few valid source pixels; skipping scene.\n');
        fail_log{end+1} = fname;
        continue
    end

    F_HH  = scatteredInterpolant(xs(src_valid),ys(src_valid),HH_raw(src_valid), ...
        'linear','none');
    F_HV  = scatteredInterpolant(xs(src_valid),ys(src_valid),HV_raw(src_valid), ...
        'linear','none');
    F_ANG = scatteredInterpolant(xs(src_valid),ys(src_valid),ANG_raw(src_valid), ...
        'linear','none');
    HH_3413  = F_HH(xq,yq);
    HV_3413  = F_HV(xq,yq);
    ANG_3413 = F_ANG(xq,yq);
    clear xs ys valid_xy src_valid F_HH F_HV F_ANG

    clear HH_raw HV_raw ANG_raw

    %% Cosine incidence-angle normalization in the linear domain
    cos_ratio  = cosd(ANG_3413) ./ cosd(theta_ref);
    HH_corr    = HH_3413 .* cos_ratio;
    HV_corr    = HV_3413 .* cos_ratio;
    HH_corr(HH_corr <= 0) = NaN;
    HV_corr(HV_corr <= 0) = NaN;

    HH_corr_dB = 10 * log10(HH_corr);
    HV_corr_dB = 10 * log10(HV_corr);

    %% Save the processed scene
    [~, stem] = fileparts(fname);
    save_path  = fullfile(output_dir, [stem '_processed.mat']);
    save(save_path, ...
         'HH_3413','HV_3413','ANG_3413', ...
         'HH_corr','HV_corr','HH_corr_dB','HV_corr_dB', ...
         'xq','yq','xv','yv','target_res','theta_ref', '-v7.3');
    fprintf('    Saved processed scene.\n');
end

%% Report failed files
fprintf('\nProcessing complete. Failed scenes: %d\n', numel(fail_log));
for i = 1:numel(fail_log)
    fprintf('    %s\n', fail_log{i});
end
