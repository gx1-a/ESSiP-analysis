%% 06_sea_ice_drift.m
% Purpose: Retrieve sea-ice displacement for 22-26 March 2026 using
% forward-backward normalized-cross-correlation template matching.
% Ridge stationarity is evaluated before ridge vectors are used to estimate
% residual pairwise coregistration offsets. Corrected displacement vectors are
% quality controlled, converted to drift velocity, and projected onto the
% ESSiP opening normal (14.39 deg true bearing).
% InternalUncertainty describes matching/registration consistency; it is not
% an error estimate relative to an independent drift observation.
% Inputs: five processed SAR scenes from 22-26 March 2026 and the rough-ice footprint file.
% MATLAB: R2022b

clear; clc; close all;


%% Repository paths
project_root = fileparts(mfilename('fullpath'));


%% ==================== 1. Input files ====================
processed_dir = fullfile(project_root,'processed');

scene_names = { ...
    'Subset_S1A_EW_GRDM_1SDH_20260322T202439_20260322T202533_063745_0803C3_25C4_Orb_NR_Cal_Spk_TC_processed.mat'; ...
    'Subset_S1C_EW_GRDM_1SDH_20260323T201534_20260323T201634_006896_00DF2D_DB8A_Orb_NR_Cal_Spk_TC_processed.mat'; ...
    'Subset_S1A_EW_GRDM_1SDH_20260324T200819_20260324T200919_063774_0804D1_7E5D_Orb_NR_Cal_Spk_TC_processed.mat'; ...
    'Subset_S1C_EW_GRDM_1SDH_20260325T195911_20260325T200011_006925_00E02C_CBC7_Orb_NR_Cal_Spk_TC_processed.mat'; ...
    'Subset_S1A_EW_GRDM_1SDH_20260326T195156_20260326T195256_063803_0805EB_3E68_Orb_NR_Cal_Spk_TC_processed.mat'};
files = cellfun(@(s) fullfile(processed_dir,s),scene_names,'UniformOutput',false);

ridge_file = fullfile(processed_dir,'ridge','edge_consensus_footprint', ...
    'ridge_edge_consensus_final.mat');
out_dir = fullfile(processed_dir,'ice_drift_event_20260322_20260326');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

proj3413 = projcrs(3413);
opening_bearing_true_deg = 14.39;  % clockwise from true north


P = struct();


P.grid_step_px = 20;              % approximately 800 m at 40 m pixels
P.template_size_px = 31;          % approximately 1.24 km x 1.24 km
P.search_radius_px = 250;         % approximately 10 km radius
P.median_filter_px = 5;           % 5 x 5 pixels


P.ncc_threshold = 0.60;
P.peak_margin_threshold = 0.05;   % primary minus secondary correlation peak
P.peak_exclusion_radius_px = 5;
P.search_edge_buffer_px = 4;
P.min_template_valid_fraction = 0.98;
P.min_target_valid_fraction = 0.98;
P.min_template_std = 0.025;
P.min_source_mean_dB = -20.0;


P.forward_backward_max_px = 3.0;


P.tracking_radius_km = 35;
P.analysis_radius_km = 30;
P.ocean_side_only = true;


P.local_filter_radius_km = 4;
P.local_min_neighbors = 4;
P.local_floor_tolerance_kmd = 0.75;
P.local_mad_multiplier = 3.0;


P.ridge_reference_min_points = 5;
P.ridge_reference_min_ncc = 0.70;
P.ridge_reference_max_fb_px = 2.5;


P.plot_stride = 2;
P.map_radius_km = 32;
P.reference_arrow_kmd = 10;


for i = 1:numel(files)
    if ~isfile(files{i})
        error('SAR file not found: %s',files{i});
    end
end

if ~isfile(ridge_file)
    error('Rough-ice-footprint file not found: %s',ridge_file);
end

R = load(ridge_file);


if isfield(R,'ridge_footprint_inner')
    ridge_mask_source = logical(R.ridge_footprint_inner);
elseif isfield(R,'ridge_footprint_consensus')
    ridge_mask_source = logical(R.ridge_footprint_consensus);
else
    error('Rough-ice file is missing ridge_footprint_inner/ridge_footprint_consensus.');
end

if ~isfield(R,'X_master') || ~isfield(R,'Y_master')
    error('Rough-ice file is missing X_master/Y_master.');
end

ridge_X = double(R.X_master);
ridge_Y = double(R.Y_master);

ridge_x = ridge_X(1,:);
ridge_y = ridge_Y(:,1);

[ridge_x,ridge_y,ridge_mask_source] = ...
    normalize_regular_grid(ridge_x,ridge_y,ridge_mask_source);

ridge_xc = mean(ridge_X(ridge_mask_source),'omitnan');
ridge_yc = mean(ridge_Y(ridge_mask_source),'omitnan');


[opening_grid_unit_x,opening_grid_unit_y] = ...
    true_bearing_to_projected_unit( ...
    ridge_xc,ridge_yc,opening_bearing_true_deg,proj3413);

[north_grid_unit_x,north_grid_unit_y] = ...
    true_bearing_to_projected_unit( ...
    ridge_xc,ridge_yc,0,proj3413);

fprintf('Rough-ice-footprint centroid: X=%.1f m, Y=%.1f m\n',ridge_xc,ridge_yc);
fprintf('True opening-normal bearing: %.2f°\n',opening_bearing_true_deg);
fprintf('14.39° unit direction in EPSG:3413: [%.4f, %.4f]\n', ...
    opening_grid_unit_x,opening_grid_unit_y);


nPair = numel(files)-1;

PairStartUTC = NaT(nPair,1,'TimeZone','UTC');
PairEndUTC = NaT(nPair,1,'TimeZone','UTC');
PairMidUTC = NaT(nPair,1,'TimeZone','UTC');
DeltaTime_h = nan(nPair,1);

NNodes = zeros(nPair,1);
NForwardAccepted = zeros(nPair,1);
NForwardBackwardAccepted = zeros(nPair,1);
NAfterLocalFilter = zeros(nPair,1);
NEventVectors = zeros(nPair,1);
NRidgeReference = zeros(nPair,1);

RegistrationDx_m = nan(nPair,1);
RegistrationDy_m = nan(nPair,1);
RegistrationMagnitude_m = nan(nPair,1);
RegistrationUncertainty_m = nan(nPair,1);
ResidualRidgeMedianSpeed_kmd = nan(nPair,1);
ResidualRidgeP95Speed_kmd = nan(nPair,1);

MedianSpeed_kmd = nan(nPair,1);
SpeedQ25_kmd = nan(nPair,1);
SpeedQ75_kmd = nan(nPair,1);

MedianOpeningComponent_kmd = nan(nPair,1);
OpeningQ25_kmd = nan(nPair,1);
OpeningQ75_kmd = nan(nPair,1);
FractionOpeningPositive = nan(nPair,1);

MedianInternalUncertainty_kmd = nan(nPair,1);
P95InternalUncertainty_kmd = nan(nPair,1);
MedianNCC = nan(nPair,1);
MedianPeakMargin = nan(nPair,1);
MedianFBError_px = nan(nPair,1);

vector_tables = cell(nPair,1);
pair_backgrounds = cell(nPair,1);
pair_xv = cell(nPair,1);
pair_yv = cell(nPair,1);
pair_ridge_masks = cell(nPair,1);


for p = 1:nPair
    fprintf('\n====================================================\n');
    fprintf('Processing image pair %d/%d\n',p,nPair);
    fprintf('%s\n%s\n',files{p},files{p+1});
    fprintf('====================================================\n');

    t1 = parse_sentinel1_midtime(files{p});
    t2 = parse_sentinel1_midtime(files{p+1});

    dt_h = hours(t2-t1);
    dt_days = dt_h/24;

    if dt_h<=0
        error('Image acquisition times are not in chronological order.');
    end

    PairStartUTC(p) = t1;
    PairEndUTC(p) = t2;
    PairMidUTC(p) = t1+(t2-t1)/2;
    DeltaTime_h(p) = dt_h;

    fprintf('Exact midpoint time difference: %.4f h\n',dt_h);


    [HH1,HH2,xv,yv] = ...
        load_pair_common_grid(files{p},files{p+1});

    [XG,YG] = meshgrid(xv,yv);

    pixel_dx = median(abs(diff(xv)),'omitnan');
    pixel_dy = median(abs(diff(yv)),'omitnan');
    pixel_res_m = mean([pixel_dx,pixel_dy]);

    fprintf('Common grid: %d x %d; pixel spacing approximately %.2f m\n', ...
        size(HH1,1),size(HH1,2),pixel_res_m);


    F_ridge = griddedInterpolant( ...
        {ridge_y,ridge_x},double(ridge_mask_source), ...
        'nearest','none');

    ridge_mask = F_ridge(YG,XG)>0.5;
    ridge_mask = imerode(ridge_mask,strel('disk',2,0));

    if ~any(ridge_mask(:))
        error('The stable rough-ice mask does not overlap the current image grid.');
    end

    pair_backgrounds{p} = HH1;
    pair_xv{p} = xv;
    pair_yv{p} = yv;
    pair_ridge_masks{p} = ridge_mask;


    [img1,valid1] = preprocess_for_matching(HH1,P);
    [img2,valid2] = preprocess_for_matching(HH2,P);


    half_t = floor(P.template_size_px/2);
    margin_px = P.search_radius_px+half_t+2;

    x_candidates = margin_px:P.grid_step_px: ...
        (size(img1,2)-margin_px);

    y_candidates = margin_px:P.grid_step_px: ...
        (size(img1,1)-margin_px);

    [XN,YN] = meshgrid(x_candidates,y_candidates);

    node_x_m = interp1(1:numel(xv),xv,XN);
    node_y_m = interp1(1:numel(yv),yv,YN);

    within_tracking = hypot( ...
        node_x_m-ridge_xc,node_y_m-ridge_yc) <= ...
        P.tracking_radius_km*1000;

    XN = XN(within_tracking);
    YN = YN(within_tracking);

    NNodes(p) = numel(XN);
    fprintf('Candidate tracking nodes: %d\n',NNodes(p));


    result = track_nodes_forward_backward( ...
        img1,img2,HH1,valid1,valid2,XN,YN,P);

    NForwardAccepted(p) = sum(result.ForwardAccepted);
    NForwardBackwardAccepted(p) = sum(result.FinalAccepted);

    fprintf('Forward QC accepted: %d\n',NForwardAccepted(p));
    fprintf('Forward-backward consistency accepted: %d\n', ...
        NForwardBackwardAccepted(p));

    if NForwardBackwardAccepted(p)<10
        warning('Fewer than 10 valid drift vectors remain for this image pair.');
    end

    keep = result.FinalAccepted;

    cx1 = result.Cx1(keep);
    cy1 = result.Cy1(keep);
    cx2 = result.Cx2(keep);
    cy2 = result.Cy2(keep);

    x1_m = interp1(1:numel(xv),xv,cx1);
    y1_m = interp1(1:numel(yv),yv,cy1);
    x2_m = interp1(1:numel(xv),xv,cx2);
    y2_m = interp1(1:numel(yv),yv,cy2);

    dx_raw_m = x2_m-x1_m;
    dy_raw_m = y2_m-y1_m;

    score = result.NCC(keep);
    peak_margin = result.PeakMargin(keep);
    fb_error_px = result.FBError_px(keep);


    in_ridge = interp2( ...
        xv,yv,double(ridge_mask), ...
        x1_m,y1_m,'nearest',0)>0.5;

    ridge_ref = in_ridge & ...
        score>=P.ridge_reference_min_ncc & ...
        fb_error_px<=P.ridge_reference_max_fb_px;

    NRidgeReference(p) = sum(ridge_ref);

    if NRidgeReference(p)>=P.ridge_reference_min_points
        reg_dx = median(dx_raw_m(ridge_ref),'omitnan');
        reg_dy = median(dy_raw_m(ridge_ref),'omitnan');

        sigma_dx = robust_sigma_of_median(dx_raw_m(ridge_ref));
        sigma_dy = robust_sigma_of_median(dy_raw_m(ridge_ref));
        reg_unc = hypot(sigma_dx,sigma_dy);
    else
        warning(['Too few high-quality reference vectors inside the rough-ice footprint; ' ...
            'no residual coregistration correction is applied to this pair.']);

        reg_dx = 0;
        reg_dy = 0;
        reg_unc = pixel_res_m/2;
    end

    RegistrationDx_m(p) = reg_dx;
    RegistrationDy_m(p) = reg_dy;
    RegistrationMagnitude_m(p) = hypot(reg_dx,reg_dy);
    RegistrationUncertainty_m(p) = reg_unc;

    dx_corr_m = dx_raw_m-reg_dx;
    dy_corr_m = dy_raw_m-reg_dy;


    U_grid_kmd = dx_corr_m/1000/dt_days;
    V_grid_kmd = dy_corr_m/1000/dt_days;

    valid_local = local_vector_filter( ...
        x1_m,y1_m,U_grid_kmd,V_grid_kmd,P);

    NAfterLocalFilter(p) = sum(valid_local);

    x1_m = x1_m(valid_local);
    y1_m = y1_m(valid_local);
    dx_raw_m = dx_raw_m(valid_local);
    dy_raw_m = dy_raw_m(valid_local);
    dx_corr_m = dx_corr_m(valid_local);
    dy_corr_m = dy_corr_m(valid_local);
    score = score(valid_local);
    peak_margin = peak_margin(valid_local);
    fb_error_px = fb_error_px(valid_local);
    in_ridge = in_ridge(valid_local);


    x2_corr_m = x1_m+dx_corr_m;
    y2_corr_m = y1_m+dy_corr_m;

    [lat1,lon1] = projinv(proj3413,x1_m,y1_m);
    [lat2,lon2] = projinv(proj3413,x2_corr_m,y2_corr_m);

    [distance_m,true_bearing_deg] = ...
        spherical_distance_bearing(lat1,lon1,lat2,lon2);

    speed_kmd = distance_m/1000/dt_days;

    U_true_east_kmd = speed_kmd.*sind(true_bearing_deg);
    V_true_north_kmd = speed_kmd.*cosd(true_bearing_deg);

    opening_component_kmd = ...
        U_true_east_kmd*sind(opening_bearing_true_deg) + ...
        V_true_north_kmd*cosd(opening_bearing_true_deg);

    cross_component_kmd = ...
        U_true_east_kmd*cosd(opening_bearing_true_deg) - ...
        V_true_north_kmd*sind(opening_bearing_true_deg);



    fb_unc_m = 0.5*fb_error_px*pixel_res_m;


    quantization_unc_m = 0.25*pixel_res_m;

    internal_uncertainty_m = sqrt( ...
        fb_unc_m.^2 + reg_unc.^2 + quantization_unc_m.^2);

    internal_uncertainty_kmd = ...
        internal_uncertainty_m/1000/dt_days;


    rel_x = x1_m-ridge_xc;
    rel_y = y1_m-ridge_yc;

    within_analysis = hypot(rel_x,rel_y) <= ...
        P.analysis_radius_km*1000;

    ocean_side_coordinate_m = ...
        rel_x*opening_grid_unit_x + ...
        rel_y*opening_grid_unit_y;

    if P.ocean_side_only
        on_event_side = ocean_side_coordinate_m>=0;
    else
        on_event_side = true(size(ocean_side_coordinate_m));
    end

    event_vector = ...
        within_analysis & on_event_side & ~in_ridge;

    NEventVectors(p) = sum(event_vector);


    if any(in_ridge)
        ResidualRidgeMedianSpeed_kmd(p) = ...
            median(speed_kmd(in_ridge),'omitnan');

        ResidualRidgeP95Speed_kmd(p) = ...
            prctile(speed_kmd(in_ridge),95);
    end


    nVec = numel(x1_m);

    PairID = repmat(p,nVec,1);
    StartUTC = repmat(t1,nVec,1);
    EndUTC = repmat(t2,nVec,1);
    Delta_h = repmat(dt_h,nVec,1);

    VectorTable = table( ...
        PairID,StartUTC,EndUTC,Delta_h, ...
        x1_m,y1_m,lat1,lon1, ...
        dx_raw_m,dy_raw_m,dx_corr_m,dy_corr_m, ...
        speed_kmd,true_bearing_deg, ...
        U_true_east_kmd,V_true_north_kmd, ...
        opening_component_kmd,cross_component_kmd, ...
        score,peak_margin,fb_error_px, ...
        internal_uncertainty_kmd, ...
        in_ridge,event_vector,ocean_side_coordinate_m, ...
        'VariableNames',{ ...
        'PairID','StartUTC','EndUTC','DeltaTime_h', ...
        'XStart_m','YStart_m','LatStart_deg','LonStart_deg', ...
        'DxRaw_m','DyRaw_m','DxCorrected_m','DyCorrected_m', ...
        'Speed_km_day','TrueBearing_deg', ...
        'UEast_km_day','VNorth_km_day', ...
        'OpeningComponent_km_day','CrossComponent_km_day', ...
        'NCC','PeakMargin','ForwardBackwardError_px', ...
        'InternalUncertainty_km_day', ...
        'InsideStableRidge','EventAnalysisVector', ...
        'OceanSideCoordinate_m'});

    vector_tables{p} = VectorTable;

    date1_str = datestr(t1,'yyyymmdd');
    date2_str = datestr(t2,'yyyymmdd');

    writetable(VectorTable,fullfile(out_dir, ...
        sprintf('SeaIceDrift_Vectors_%s_%s.csv', ...
        date1_str,date2_str)));

    save(fullfile(out_dir, ...
        sprintf('SeaIceDrift_Vectors_%s_%s.mat', ...
        date1_str,date2_str)), ...
        'VectorTable','HH1','HH2','xv','yv','ridge_mask','P');


    if any(event_vector)
        MedianSpeed_kmd(p) = ...
            median(speed_kmd(event_vector),'omitnan');

        SpeedQ25_kmd(p) = ...
            prctile(speed_kmd(event_vector),25);

        SpeedQ75_kmd(p) = ...
            prctile(speed_kmd(event_vector),75);

        MedianOpeningComponent_kmd(p) = ...
            median(opening_component_kmd(event_vector),'omitnan');

        OpeningQ25_kmd(p) = ...
            prctile(opening_component_kmd(event_vector),25);

        OpeningQ75_kmd(p) = ...
            prctile(opening_component_kmd(event_vector),75);

        FractionOpeningPositive(p) = ...
            mean(opening_component_kmd(event_vector)>0);

        MedianInternalUncertainty_kmd(p) = ...
            median(internal_uncertainty_kmd(event_vector),'omitnan');

        P95InternalUncertainty_kmd(p) = ...
            prctile(internal_uncertainty_kmd(event_vector),95);

        MedianNCC(p) = median(score(event_vector),'omitnan');
        MedianPeakMargin(p) = ...
            median(peak_margin(event_vector),'omitnan');
        MedianFBError_px(p) = ...
            median(fb_error_px(event_vector),'omitnan');
    end

    fprintf('Stable-footprint registration reference points: %d\n',NRidgeReference(p));
    fprintf('Residual coregistration offset: dx=%.1f m, dy=%.1f m, |d|=%.1f m\n', ...
        reg_dx,reg_dy,hypot(reg_dx,reg_dy));
    fprintf('Vectors after local consistency filtering: %d\n',NAfterLocalFilter(p));
    fprintf('Opening-side event vectors: %d\n',NEventVectors(p));
    fprintf('Median drift speed: %.2f km d^-1\n',MedianSpeed_kmd(p));
    fprintf('Median component along 14.39°: %.2f km d^-1\n', ...
        MedianOpeningComponent_kmd(p));
    fprintf('Fraction moving in the opening direction: %.1f%%\n', ...
        100*FractionOpeningPositive(p));
    fprintf('Median internal uncertainty: %.3f km d^-1\n', ...
        MedianInternalUncertainty_kmd(p));
end


PairSummary = table( ...
    (1:nPair)',PairStartUTC,PairEndUTC,PairMidUTC,DeltaTime_h, ...
    NNodes,NForwardAccepted,NForwardBackwardAccepted, ...
    NAfterLocalFilter,NEventVectors,NRidgeReference, ...
    RegistrationDx_m,RegistrationDy_m, ...
    RegistrationMagnitude_m,RegistrationUncertainty_m, ...
    ResidualRidgeMedianSpeed_kmd,ResidualRidgeP95Speed_kmd, ...
    MedianSpeed_kmd,SpeedQ25_kmd,SpeedQ75_kmd, ...
    MedianOpeningComponent_kmd,OpeningQ25_kmd,OpeningQ75_kmd, ...
    FractionOpeningPositive, ...
    MedianInternalUncertainty_kmd,P95InternalUncertainty_kmd, ...
    MedianNCC,MedianPeakMargin,MedianFBError_px, ...
    'VariableNames',{ ...
    'PairID','StartUTC','EndUTC','MidUTC','DeltaTime_h', ...
    'NNodes','NForwardAccepted','NForwardBackwardAccepted', ...
    'NAfterLocalFilter','NEventVectors','NRidgeReference', ...
    'RegistrationDx_m','RegistrationDy_m', ...
    'RegistrationMagnitude_m','RegistrationUncertainty_m', ...
    'ResidualRidgeMedianSpeed_km_day', ...
    'ResidualRidgeP95Speed_km_day', ...
    'MedianSpeed_km_day','SpeedQ25_km_day','SpeedQ75_km_day', ...
    'MedianOpeningComponent_km_day', ...
    'OpeningQ25_km_day','OpeningQ75_km_day', ...
    'FractionOpeningPositive', ...
    'MedianInternalUncertainty_km_day', ...
    'P95InternalUncertainty_km_day', ...
    'MedianNCC','MedianPeakMargin', ...
    'MedianForwardBackwardError_px'});

fprintf('\n===== Consecutive-event sea-ice drift summary =====\n');
disp(PairSummary);

writetable(PairSummary,fullfile(out_dir, ...
    'SeaIceDrift_Event_PairSummary.csv'));

save(fullfile(out_dir,'SeaIceDrift_Event_AllResults.mat'), ...
    'PairSummary','vector_tables','P', ...
    'opening_bearing_true_deg','scene_names','-v7.3');

%% ==================== 7. Finish ====================
fprintf('\nSea-ice drift outputs saved to: %s\n',out_dir);



function mid_time = parse_sentinel1_midtime(file_path)
    [~,name,~] = fileparts(file_path);

    tokens = regexp(name, ...
        '(\d{8}T\d{6})_(\d{8}T\d{6})', ...
        'tokens','once');

    if isempty(tokens)
        error('Could not parse acquisition time from filename: %s',name);
    end

    start_time = datetime(tokens{1}, ...
        'InputFormat','yyyyMMdd''T''HHmmss', ...
        'TimeZone','UTC');

    end_time = datetime(tokens{2}, ...
        'InputFormat','yyyyMMdd''T''HHmmss', ...
        'TimeZone','UTC');

    mid_time = start_time+(end_time-start_time)/2;
end

function [HH1,HH2,xv,yv] = load_pair_common_grid(file1,file2)
    A = load(file1,'HH_corr_dB','xv','yv');
    B = load(file2,'HH_corr_dB','xv','yv');

    [x1,y1,HH1] = normalize_regular_grid( ...
        double(A.xv(:)'),double(A.yv(:)), ...
        double(A.HH_corr_dB));

    [x2,y2,HH2_original] = normalize_regular_grid( ...
        double(B.xv(:)'),double(B.yv(:)), ...
        double(B.HH_corr_dB));

    xv = x1;
    yv = y1;

    same_x = numel(x1)==numel(x2) && ...
        max(abs(x1-x2),[],'omitnan')<1e-6;

    same_y = numel(y1)==numel(y2) && ...
        max(abs(y1-y2),[],'omitnan')<1e-6;

    if same_x && same_y
        HH2 = HH2_original;
    else
        fprintf('Image grids differ; interpolating the second scene to the first-scene grid.\n');

        F = griddedInterpolant( ...
            {y2,x2},HH2_original,'linear','none');

        [X,Y] = meshgrid(xv,yv);
        HH2 = F(Y,X);
    end
end

function [x,y,A] = normalize_regular_grid(x,y,A)
    x = double(x(:)');
    y = double(y(:));

    if numel(x)~=size(A,2) || numel(y)~=size(A,1)
        error('Coordinate lengths do not match image dimensions.');
    end

    if x(1)>x(end)
        x = fliplr(x);
        A = fliplr(A);
    end

    if y(1)>y(end)
        y = flipud(y);
        A = flipud(A);
    end
end

function [img,valid] = preprocess_for_matching(HH,P)
    valid = isfinite(HH);

    img = mat2gray(HH,[-25 0]);
    img(~valid) = 0;

    img = medfilt2(img, ...
        [P.median_filter_px P.median_filter_px], ...
        'symmetric');
end

function result = track_nodes_forward_backward( ...
    img1,img2,HH1,valid1,valid2,XN,YN,P)

    n = numel(XN);

    Cx1 = double(XN(:));
    Cy1 = double(YN(:));
    Cx2 = nan(n,1);
    Cy2 = nan(n,1);
    NCC = nan(n,1);
    PeakMargin = nan(n,1);
    FBError_px = nan(n,1);

    ForwardAccepted = false(n,1);
    FinalAccepted = false(n,1);

    fprintf('Running forward-backward NCC matching for %d nodes...\n',n);

    for i = 1:n
        cx = Cx1(i);
        cy = Cy1(i);

        [ok_f,cx2,cy2,score,margin] = ...
            match_one_point( ...
            img1,img2,HH1,valid1,valid2,cx,cy,P);

        if ~ok_f
            continue;
        end

        ForwardAccepted(i) = true;
        Cx2(i) = cx2;
        Cy2(i) = cy2;
        NCC(i) = score;
        PeakMargin(i) = margin;

        cx2_int = round(cx2);
        cy2_int = round(cy2);

        [ok_b,cx_back,cy_back,~,~] = ...
            match_one_point( ...
            img2,img1,[],valid2,valid1, ...
            cx2_int,cy2_int,P);

        if ~ok_b
            continue;
        end

        closure = hypot(cx_back-cx,cy_back-cy);
        FBError_px(i) = closure;

        if closure<=P.forward_backward_max_px
            FinalAccepted(i) = true;
        end
    end

    result = table( ...
        Cx1,Cy1,Cx2,Cy2,NCC,PeakMargin,FBError_px, ...
        ForwardAccepted,FinalAccepted);
end

function [ok,match_cx,match_cy,max_score,peak_margin] = ...
    match_one_point( ...
    imgA,imgB,HHA,validA,validB,cx,cy,P)

    ok = false;
    match_cx = NaN;
    match_cy = NaN;
    max_score = NaN;
    peak_margin = NaN;

    half_t = floor(P.template_size_px/2);
    R = P.search_radius_px;

    if cx-half_t<1 || cx+half_t>size(imgA,2) || ...
       cy-half_t<1 || cy+half_t>size(imgA,1)
        return;
    end

    if cx-R<1 || cx+R>size(imgB,2) || ...
       cy-R<1 || cy+R>size(imgB,1)
        return;
    end

    template_valid = validA( ...
        cy-half_t:cy+half_t, ...
        cx-half_t:cx+half_t);

    if mean(template_valid(:))<P.min_template_valid_fraction
        return;
    end

    if ~isempty(HHA)
        source_patch_dB = HHA( ...
            cy-half_t:cy+half_t, ...
            cx-half_t:cx+half_t);

        if mean(source_patch_dB,'all','omitnan') < ...
                P.min_source_mean_dB
            return;
        end
    end

    template = imgA( ...
        cy-half_t:cy+half_t, ...
        cx-half_t:cx+half_t);

    if std(template,0,'all','omitnan')<P.min_template_std
        return;
    end

    search_region = imgB( ...
        cy-R:cy+R, ...
        cx-R:cx+R);

    search_valid = validB( ...
        cy-R:cy+R, ...
        cx-R:cx+R);

    ncc = normxcorr2(template,search_region);

    t = P.template_size_px;


    valid_ncc = ncc(t:end-t+1,t:end-t+1);

    kernel = ones(t,t);

    target_valid_fraction = ...
        conv2(double(search_valid),kernel,'valid')/(t*t);

    valid_ncc(target_valid_fraction < ...
        P.min_target_valid_fraction) = -Inf;

    [max_score,max_idx] = max(valid_ncc(:));

    if ~isfinite(max_score) || max_score<P.ncc_threshold
        return;
    end

    [py,px] = ind2sub(size(valid_ncc),max_idx);

    if px<=P.search_edge_buffer_px || ...
       py<=P.search_edge_buffer_px || ...
       px>size(valid_ncc,2)-P.search_edge_buffer_px || ...
       py>size(valid_ncc,1)-P.search_edge_buffer_px
        return;
    end

    secondary = valid_ncc;

    rr1 = max(1,py-P.peak_exclusion_radius_px);
    rr2 = min(size(secondary,1), ...
        py+P.peak_exclusion_radius_px);

    cc1 = max(1,px-P.peak_exclusion_radius_px);
    cc2 = min(size(secondary,2), ...
        px+P.peak_exclusion_radius_px);

    secondary(rr1:rr2,cc1:cc2) = -Inf;
    second_score = max(secondary(:));

    if isfinite(second_score)
        peak_margin = max_score-second_score;

        if peak_margin<P.peak_margin_threshold
            return;
        end
    else
        peak_margin = Inf;
    end

    [sub_dx,sub_dy] = ...
        parabolic_subpixel_peak(valid_ncc,px,py);



    match_cx = cx-R+(px-1)+half_t+sub_dx;
    match_cy = cy-R+(py-1)+half_t+sub_dy;

    ok = true;
end

function [sub_dx,sub_dy] = ...
    parabolic_subpixel_peak(A,px,py)

    sub_dx = 0;
    sub_dy = 0;

    if px>1 && px<size(A,2)
        left = A(py,px-1);
        center = A(py,px);
        right = A(py,px+1);

        denom = left-2*center+right;

        if isfinite(denom) && abs(denom)>eps
            sub_dx = 0.5*(left-right)/denom;
            sub_dx = max(-1,min(1,sub_dx));
        end
    end

    if py>1 && py<size(A,1)
        up = A(py-1,px);
        center = A(py,px);
        down = A(py+1,px);

        denom = up-2*center+down;

        if isfinite(denom) && abs(denom)>eps
            sub_dy = 0.5*(up-down)/denom;
            sub_dy = max(-1,min(1,sub_dy));
        end
    end
end

function keep = local_vector_filter(x,y,U,V,P)
    n = numel(x);
    keep = true(n,1);

    radius_m = P.local_filter_radius_km*1000;

    for i = 1:n
        d = hypot(x-x(i),y-y(i));
        neighbors = d<=radius_m;

        if sum(neighbors)<P.local_min_neighbors
            continue;
        end

        medU = median(U(neighbors),'omitnan');
        medV = median(V(neighbors),'omitnan');

        residual_neighbors = hypot( ...
            U(neighbors)-medU, ...
            V(neighbors)-medV);

        med_residual = median(residual_neighbors,'omitnan');

        robust_scale = 1.4826*median( ...
            abs(residual_neighbors-med_residual), ...
            'omitnan');

        tolerance = max( ...
            P.local_floor_tolerance_kmd, ...
            P.local_mad_multiplier*robust_scale);

        residual_i = hypot(U(i)-medU,V(i)-medV);

        if residual_i>tolerance
            keep(i) = false;
        end
    end
end

function sigma_median = robust_sigma_of_median(x)
    x = x(isfinite(x));

    if isempty(x)
        sigma_median = NaN;
        return;
    end

    medx = median(x);
    robust_sigma = 1.4826*median(abs(x-medx));

    sigma_median = ...
        1.2533*robust_sigma/sqrt(numel(x));

    if ~isfinite(sigma_median)
        sigma_median = 0;
    end
end

function [distance_m,bearing_deg] = ...
    spherical_distance_bearing(lat1,lon1,lat2,lon2)

    earth_radius_m = 6371008.8;

    phi1 = deg2rad(lat1);
    phi2 = deg2rad(lat2);
    dphi = phi2-phi1;
    dlambda = deg2rad(lon2-lon1);

    a = sin(dphi/2).^2 + ...
        cos(phi1).*cos(phi2).*sin(dlambda/2).^2;

    c = 2*atan2(sqrt(a),sqrt(max(0,1-a)));
    distance_m = earth_radius_m*c;

    y = sin(dlambda).*cos(phi2);
    x = cos(phi1).*sin(phi2) - ...
        sin(phi1).*cos(phi2).*cos(dlambda);

    bearing_deg = mod(atan2d(y,x),360);
    bearing_deg(distance_m<1e-9) = 0;
end

function [ux,uy] = true_bearing_to_projected_unit( ...
    x0,y0,bearing_deg,projcrs_object)

    [lat0,lon0] = projinv(projcrs_object,x0,y0);

    distance_m = 1000;
    earth_radius_m = 6371008.8;
    delta = distance_m/earth_radius_m;
    theta = deg2rad(bearing_deg);

    phi1 = deg2rad(lat0);
    lambda1 = deg2rad(lon0);

    phi2 = asin( ...
        sin(phi1)*cos(delta) + ...
        cos(phi1)*sin(delta)*cos(theta));

    lambda2 = lambda1+atan2( ...
        sin(theta)*sin(delta)*cos(phi1), ...
        cos(delta)-sin(phi1)*sin(phi2));

    lat2 = rad2deg(phi2);
    lon2 = rad2deg(lambda2);

    [x2,y2] = projfwd(projcrs_object,lat2,lon2);

    ux = x2-x0;
    uy = y2-y0;

    scale = hypot(ux,uy);
    ux = ux/scale;
    uy = uy/scale;
end
function [HH_geo, ridge_geo, lon_vec, lat_vec] = ...
    project_pair_to_geographic_grid(HH, ridge_mask, xv, yv, proj3413, nGrid)

    [X, Y] = meshgrid(xv, yv);
    [LAT, LON] = projinv(proj3413, X, Y);

    valid = isfinite(HH) & isfinite(LAT) & isfinite(LON);

    lon_min = min(LON(valid));
    lon_max = max(LON(valid));
    lat_min = min(LAT(valid));
    lat_max = max(LAT(valid));

    lon_vec = linspace(lon_min, lon_max, nGrid);
    lat_vec = linspace(lat_min, lat_max, nGrid);
    [LONq, LATq] = meshgrid(lon_vec, lat_vec);

    F_hh = scatteredInterpolant(LON(valid), LAT(valid), HH(valid), ...
        'linear', 'none');
    HH_geo = F_hh(LONq, LATq);

    F_ridge = scatteredInterpolant(LON(valid), LAT(valid), ...
        double(ridge_mask(valid)), 'nearest', 'none');
    ridge_geo = F_ridge(LONq, LATq) > 0.5;
end
function [lat2, lon2] = move_point_geographic(lat1, lon1, bearing_deg, distance_km)
    R = 6371.0088; % Earth radius, km

    phi1 = deg2rad(lat1);
    lambda1 = deg2rad(lon1);
    theta = deg2rad(bearing_deg);
    delta = distance_km / R;

    phi2 = asin( sin(phi1).*cos(delta) + cos(phi1).*sin(delta).*cos(theta) );
    lambda2 = lambda1 + atan2( sin(theta).*sin(delta).*cos(phi1), ...
                               cos(delta) - sin(phi1).*sin(phi2) );

    lat2 = rad2deg(phi2);
    lon2 = rad2deg(lambda2);
end
