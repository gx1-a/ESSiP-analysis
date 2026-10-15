%% 01_derive_rough_ice_footprint.m
% Derive the persistent rough-ice footprint used in the manuscript.
%
% This script implements the complete manuscript/SI workflow in one file:
%   1) calculate multitemporal brightness-anomaly and texture persistence
%      from the 100 Sentinel-1 HH scenes;
%   2) identify stable high-persistence core candidates;
%   3) group adjacent core fragments and select the highest-scoring group;
%   4) delineate the full footprint with an edge-constrained active-contour
%      ensemble; and
%   5) retain the 50% consensus footprint as the final manuscript footprint.
%
% The 75% inner consensus is also saved for the sea-ice-motion workflow.
% The 25% outer consensus is a diagnostic only.
%
% MATLAB R2022b. Requires Image Processing Toolbox and Statistics and Machine Learning Toolbox.

clear; clc; close all;

%% 1. Repository paths and common settings
project_root = fileparts(mfilename('fullpath'));
processed_dir = fullfile(project_root,'processed');
time_csv = fullfile(processed_dir, ...
    'Sentinel1_acquisition_times_20251122_20260430.csv');
out_dir = fullfile(processed_dir,'ridge','edge_consensus_footprint');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

pixel_res = 40;                 % m
pixel_area_km2 = pixel_res^2/1e6;

%% 2. Per-scene feature parameters
small_smooth_sigma_px = 1;      % suppress speckle-scale variability
background_sigma_px = 25;       % approximately 1 km local background
texture_window = 5;             % 5 x 5 local standard deviation


bright_z_threshold = 0.8;
texture_z_threshold = 0.8;


invalid_edge_buffer_px = 3;


statistics_sample_stride = 4;

%% 3. Stable-core detection parameters





core_score_percentile_list = [99.9 99.7 99.5 99.0 98.5];


core_score_weights = [0.30 0.40 0.15 0.15];



min_core_area_km2 = 0.05;
core_close_radius_px = 2;
core_open_radius_px = 0;

%% 4. Read the 100-scene acquisition manifest
if ~isfile(time_csv)
    error('Acquisition manifest not found: %s',time_csv);
end

T = readtable(time_csv,'TextType','string');

if ~ismember('Name',T.Properties.VariableNames)
    error('Acquisition manifest is missing the Name column.');
end

file_names = string(T.Name);

if numel(file_names) ~= 100
    warning('Acquisition manifest contains %d scenes; 100 were expected.',numel(file_names));
end

if numel(unique(file_names)) ~= numel(file_names)
    error('The acquisition manifest contains duplicate filenames.');
end

fprintf('SAR scenes in acquisition manifest: %d\n',numel(file_names));

%% 5. Establish the common grid and accumulators
first_file = fullfile(processed_dir,file_names(1));
G = load(first_file,'xq','yq');

X_master = G.xq;
Y_master = G.yq;
[rows,cols] = size(X_master);

valid_count = zeros(rows,cols,'uint16');
bright_rough_count = zeros(rows,cols,'uint16');
rough_count = zeros(rows,cols,'uint16');

sum_anomaly_z = zeros(rows,cols,'single');
sum_texture_z = zeros(rows,cols,'single');
sum_HH = zeros(rows,cols,'single');

fprintf('Common grid: %d x %d\n',rows,cols);

%% 6. Calculate and accumulate scene-wise features
for i = 1:numel(file_names)
    source_file = fullfile(processed_dir,file_names(i));

    fprintf('[%3d/%d] %s\n',i,numel(file_names),file_names(i));

    if ~isfile(source_file)
        error('Source file listed in the manifest does not exist: %s',source_file);
    end

    S = load(source_file,'xq','yq','HH_corr_dB');

    if ~isfield(S,'HH_corr_dB') || ~isfield(S,'xq') || ~isfield(S,'yq')
        error('%s is missing HH_corr_dB/xq/yq.',source_file);
    end

    HH = single(S.HH_corr_dB);


    if isequal(size(S.xq),size(X_master))
        coord_diff = ...
            max(abs(double(S.xq(:))-double(X_master(:)))) + ...
            max(abs(double(S.yq(:))-double(Y_master(:))));
    else
        coord_diff = inf;
    end

    if isequal(size(HH),[rows cols]) && coord_diff < 20
        HH_aligned = HH;
    else
        HH_aligned = single(interp2( ...
            double(S.xq),double(S.yq),double(HH), ...
            double(X_master),double(Y_master), ...
            'linear',NaN));
    end

    valid_raw = isfinite(HH_aligned);

    if nnz(valid_raw)<1000
        warning('%s has too few valid pixels and was skipped.',file_names(i));
        continue;
    end


    HH_clean = HH_aligned;
    fill_value = median(HH_aligned(valid_raw),'omitnan');

    if ~isfinite(fill_value)
        fill_value = -15;
    end

    HH_clean(~valid_raw) = fill_value;
    HH_clean(HH_clean < -35) = -35;
    HH_clean(HH_clean > 0) = 0;

    HH_small = imgaussfilt(HH_clean,small_smooth_sigma_px, ...
        'Padding','replicate');

    HH_background = imgaussfilt(HH_small,background_sigma_px, ...
        'Padding','replicate');

    local_anomaly = HH_small-HH_background;
    local_texture = stdfilt(HH_small,ones(texture_window));

    valid_support = imerode(valid_raw, ...
        strel('square',2*invalid_edge_buffer_px+1));


    valid_idx = find(valid_support);
    sample_idx = valid_idx(1:statistics_sample_stride:end);

    anomaly_sample = double(local_anomaly(sample_idx));
    texture_sample = double(local_texture(sample_idx));

    [anomaly_center,anomaly_scale] = robust_center_scale(anomaly_sample);
    [texture_center,texture_scale] = robust_center_scale(texture_sample);

    anomaly_z = (double(local_anomaly)-anomaly_center)/anomaly_scale;
    texture_z = (double(local_texture)-texture_center)/texture_scale;

    bright_rough_hit = valid_support & ...
        anomaly_z>bright_z_threshold & ...
        texture_z>texture_z_threshold;

    rough_hit = valid_support & ...
        texture_z>texture_z_threshold;

    valid_count(valid_support) = valid_count(valid_support)+1;
    bright_rough_count(bright_rough_hit) = ...
        bright_rough_count(bright_rough_hit)+1;
    rough_count(rough_hit) = rough_count(rough_hit)+1;

    sum_anomaly_z(valid_support) = ...
        sum_anomaly_z(valid_support)+single(anomaly_z(valid_support));

    sum_texture_z(valid_support) = ...
        sum_texture_z(valid_support)+single(texture_z(valid_support));

    sum_HH(valid_support) = ...
        sum_HH(valid_support)+HH_aligned(valid_support);

    clear S HH HH_aligned HH_clean HH_small HH_background
    clear local_anomaly local_texture anomaly_z texture_z
end

%% 7. Build multitemporal persistence fields
valid_count_double = double(valid_count);
has_obs = valid_count_double>0;

P_bright_rough = nan(rows,cols);
P_rough = nan(rows,cols);
mean_anomaly_z = nan(rows,cols);
mean_texture_z = nan(rows,cols);
mean_HH = nan(rows,cols);

P_bright_rough(has_obs) = ...
    double(bright_rough_count(has_obs))./valid_count_double(has_obs);

P_rough(has_obs) = ...
    double(rough_count(has_obs))./valid_count_double(has_obs);

mean_anomaly_z(has_obs) = ...
    double(sum_anomaly_z(has_obs))./valid_count_double(has_obs);

mean_texture_z(has_obs) = ...
    double(sum_texture_z(has_obs))./valid_count_double(has_obs);

mean_HH(has_obs) = ...
    double(sum_HH(has_obs))./valid_count_double(has_obs);


enough_obs = valid_count_double >= 0.90*numel(file_names);

fprintf('Area with valid observations in >=90%% of scenes: %.2f%%\n', ...
    100*nnz(enough_obs)/numel(enough_obs));

%% 8. Detect stable core candidates

p_list = [90 95 97.5 99 99.5 99.9 100];

fprintf('\n===== Persistence diagnostics =====\n');
fprintf('P_bright_rough percentiles:\n');
disp(array2table(prctile(P_bright_rough(enough_obs),p_list), ...
    'VariableNames',compose('P_%g',p_list)));

fprintf('P_rough percentiles:\n');
disp(array2table(prctile(P_rough(enough_obs),p_list), ...
    'VariableNames',compose('P_%g',p_list)));

fprintf('mean_anomaly_z percentiles:\n');
disp(array2table(prctile(mean_anomaly_z(enough_obs),p_list), ...
    'VariableNames',compose('P_%g',p_list)));

fprintf('mean_texture_z percentiles:\n');
disp(array2table(prctile(mean_texture_z(enough_obs),p_list), ...
    'VariableNames',compose('P_%g',p_list)));


anomaly_support = 1./(1+exp(-mean_anomaly_z));
texture_support = 1./(1+exp(-mean_texture_z));

w = core_score_weights/sum(core_score_weights);

core_score = ...
    w(1)*P_bright_rough + ...
    w(2)*P_rough + ...
    w(3)*anomaly_support + ...
    w(4)*texture_support;

core_score(~enough_obs) = NaN;

fprintf('core_score range: %.4f–%.4f\n', ...
    min(core_score(enough_obs),[],'omitnan'), ...
    max(core_score(enough_obs),[],'omitnan'));

min_core_pixels = ceil(min_core_area_km2/pixel_area_km2);

core_candidates = false(rows,cols);
chosen_core_percentile = NaN;
chosen_core_threshold = NaN;

for pp = core_score_percentile_list
    threshold_pp = prctile(core_score(enough_obs),pp);

    test_core = enough_obs & core_score>=threshold_pp;

    if core_close_radius_px>0
        test_core = imclose(test_core, ...
            strel('disk',core_close_radius_px,0));
    end

    if core_open_radius_px>0
        test_core = imopen(test_core, ...
            strel('disk',core_open_radius_px,0));
    end

    n_before_area_filter = nnz(test_core);

    test_candidates = bwareafilt( ...
        test_core,[min_core_pixels inf],8);

    CC_test = bwconncomp(test_candidates,8);

    fprintf(['Core percentile %.1f%%: threshold=%.4f, ' ...
        '%d pixels before area filtering, %d candidate objects\n'], ...
        pp,threshold_pp,n_before_area_filter,CC_test.NumObjects);

    if CC_test.NumObjects>0
        core_candidates = test_candidates;
        chosen_core_percentile = pp;
        chosen_core_threshold = threshold_pp;
        break;
    end
end

CC = bwconncomp(core_candidates,8);

if CC.NumObjects==0
    error(['No core candidates remain after testing down to the %.1f percentile. ' ...
        'This usually indicates that min_core_area_km2 is too large or the ' ...
        'scene-wise brightness/texture thresholds are too strict. Consider reducing ' ...
        'min_core_area_km2 from %.3f to 0.01 or lowering bright_z_threshold and ' ...
        'texture_z_threshold from %.2f/%.2f to 0.5.'], ...
        core_score_percentile_list(end),min_core_area_km2, ...
        bright_z_threshold,texture_z_threshold);
end

fprintf('Selected core percentile: %.1f%%; threshold: %.4f\n', ...
    chosen_core_percentile,chosen_core_threshold);

stats = regionprops(CC, ...
    'Area','Centroid','Solidity','Perimeter','BoundingBox');

CandidateID = (1:CC.NumObjects)';
CoreArea_km2 = nan(CC.NumObjects,1);
MeanCoreScore = nan(CC.NumObjects,1);
PeakCoreScore = nan(CC.NumObjects,1);
MeanBrightRoughPersistence = nan(CC.NumObjects,1);
MeanRoughPersistence = nan(CC.NumObjects,1);
MeanAnomalyZ = nan(CC.NumObjects,1);
MeanTextureZ = nan(CC.NumObjects,1);
PersistentMass_km2 = nan(CC.NumObjects,1);
Solidity = nan(CC.NumObjects,1);
CentroidX_m = nan(CC.NumObjects,1);
CentroidY_m = nan(CC.NumObjects,1);
SelectionScore = nan(CC.NumObjects,1);

for k = 1:CC.NumObjects
    pix = CC.PixelIdxList{k};

    CoreArea_km2(k) = numel(pix)*pixel_area_km2;
    MeanCoreScore(k) = mean(core_score(pix),'omitnan');
    PeakCoreScore(k) = max(core_score(pix),[],'omitnan');

    MeanBrightRoughPersistence(k) = ...
        mean(P_bright_rough(pix),'omitnan');

    MeanRoughPersistence(k) = ...
        mean(P_rough(pix),'omitnan');

    MeanAnomalyZ(k) = ...
        mean(mean_anomaly_z(pix),'omitnan');

    MeanTextureZ(k) = ...
        mean(mean_texture_z(pix),'omitnan');


    PersistentMass_km2(k) = ...
        sum(core_score(pix),'omitnan')*pixel_area_km2;

    Solidity(k) = stats(k).Solidity;

    ccol = max(1,min(cols,round(stats(k).Centroid(1))));
    crow = max(1,min(rows,round(stats(k).Centroid(2))));

    CentroidX_m(k) = double(X_master(crow,ccol));
    CentroidY_m(k) = double(Y_master(crow,ccol));

    SelectionScore(k) = ...
        PersistentMass_km2(k) * ...
        (0.7+0.3*Solidity(k));
end

CandidateTable = table( ...
    CandidateID,CoreArea_km2,MeanCoreScore,PeakCoreScore, ...
    MeanBrightRoughPersistence,MeanRoughPersistence, ...
    MeanAnomalyZ,MeanTextureZ,PersistentMass_km2, ...
    Solidity,CentroidX_m,CentroidY_m,SelectionScore);

CandidateTable = sortrows(CandidateTable, ...
    'SelectionScore','descend');

fprintf('\n===== Persistent core candidates =====\n');
disp(CandidateTable);


%% 9. Edge-constrained footprint parameters
% Group adjacent stable-core fragments after 1.2-km dilation.
core_link_radius_m = 1200;

% Local processing region: convex hull of the selected core group + 6 km.
crop_margin_km = 6.0;

% Edge-control field: [mean brightness anomaly, mean texture, roughness persistence].
feature_weights = [0.50 0.30 0.20];
feature_smooth_sigma_px = 2;

% Nine active-contour parameter combinations.
initial_dilation_km_list = [0.3 0.6 0.9];
contraction_bias_list = [-0.25 -0.15 -0.05];
n_iterations = 350;
smooth_factor = 1.0;

% Postprocessing.
close_radius_px = 8;
open_radius_px = 3;
min_final_area_km2 = 1.0;

% Consensus thresholds.
consensus_main_fraction = 0.50;
consensus_inner_fraction = 0.75;
consensus_outer_fraction = 0.25;

% Package persistence fields for the footprint stage.
S = struct();
S.core_score = core_score;
S.P_rough = P_rough;
S.mean_anomaly_z = mean_anomaly_z;
S.mean_texture_z = mean_texture_z;
S.mean_HH = mean_HH;
S.valid_count = valid_count;
S.X_master = X_master;
S.Y_master = Y_master;
S.file_names = file_names;
S.parameters = struct('chosen_core_threshold',chosen_core_threshold, ...
    'chosen_core_percentile',chosen_core_percentile);

%% 10. Rebuild and group stable core fragments
core_candidates = enough_obs & ...
    S.core_score>=chosen_core_threshold;

core_candidates = imclose(core_candidates, ...
    strel('disk',2,0));

min_core_pixels = ceil(0.05/pixel_area_km2);

core_candidates = bwareafilt( ...
    core_candidates,[min_core_pixels inf],8);

if ~any(core_candidates(:))
    error('Stable core candidates could not be reconstructed.');
end

link_radius_px = round(core_link_radius_m/pixel_res);

linked = imdilate(core_candidates, ...
    strel('disk',link_radius_px,0));

CC_group = bwconncomp(linked,8);

group_masks = cell(CC_group.NumObjects,1);
group_scores = zeros(CC_group.NumObjects,1);
group_fragments = zeros(CC_group.NumObjects,1);
group_area_km2 = zeros(CC_group.NumObjects,1);

for g = 1:CC_group.NumObjects
    zone = false(rows,cols);
    zone(CC_group.PixelIdxList{g}) = true;

    original_core = core_candidates & zone;
    group_masks{g} = original_core;

    CC_frag = bwconncomp(original_core,8);
    group_fragments(g) = CC_frag.NumObjects;

    pix = find(original_core);
    group_area_km2(g) = numel(pix)*pixel_area_km2;

    score_mass = ...
        sum(S.core_score(pix),'omitnan')*pixel_area_km2;

    group_scores(g) = ...
        score_mass*(1+0.05*log1p(group_fragments(g)));
end

[~,selected_group_id] = max(group_scores);
seed_union = group_masks{selected_group_id};
seed_hull = bwconvhull(seed_union);

GroupID = (1:CC_group.NumObjects)';
CoreArea_km2_Group = group_area_km2;
FragmentCount = group_fragments;
GroupScore = group_scores;
Selected = GroupID==selected_group_id;
CoreGroupTable = table(GroupID,CoreArea_km2_Group,FragmentCount,GroupScore,Selected);

fprintf('Selected core group: %d; fragments: %d; core area: %.3f km2\n', ...
    selected_group_id, ...
    group_fragments(selected_group_id), ...
    group_area_km2(selected_group_id));

%% 11. Define the local processing region
crop_margin_px = round(crop_margin_km*1000/pixel_res);

crop_roi = imdilate(seed_hull, ...
    strel('disk',crop_margin_px,0));

crop_roi = crop_roi & enough_obs;

[r_idx,c_idx] = find(crop_roi);

r1 = max(1,min(r_idx));
r2 = min(rows,max(r_idx));
c1 = max(1,min(c_idx));
c2 = min(cols,max(c_idx));

crop_roi_c = crop_roi(r1:r2,c1:c2);
seed_union_c = seed_union(r1:r2,c1:c2);
seed_hull_c = seed_hull(r1:r2,c1:c2);

%% 12. Build the edge-control field
A = robust_rescale01(S.mean_anomaly_z,crop_roi,1,99.5);
T = robust_rescale01(S.mean_texture_z,crop_roi,1,99.5);
R = robust_rescale01(S.P_rough,crop_roi,1,99.5);

w = feature_weights/sum(feature_weights);

edge_feature = ...
    w(1)*A+w(2)*T+w(3)*R;

edge_feature = imgaussfilt(edge_feature, ...
    feature_smooth_sigma_px,'Padding','replicate');

edge_feature(~crop_roi) = 0;

edge_feature_c = edge_feature(r1:r2,c1:c2);


valid_values = edge_feature_c(crop_roi_c);
lo = prctile(valid_values,1);
hi = prctile(valid_values,99.5);

edge_feature_c = (edge_feature_c-lo)/max(hi-lo,eps);
edge_feature_c = max(0,min(1,edge_feature_c));
edge_feature_c(~crop_roi_c) = 0;

%% 13. Edge-constrained active-contour ensemble
nCase = numel(initial_dilation_km_list)* ...
        numel(contraction_bias_list);

InitialDilation_km = nan(nCase,1);
ContractionBias = nan(nCase,1);
Area_km2 = nan(nCase,1);
MajorAxisEllipse_km = nan(nCase,1);
MinorAxisEllipse_km = nan(nCase,1);
AlongAxisSpan_km = nan(nCase,1);
AxisAzimuth_deg = nan(nCase,1);
TouchesCropBoundary = false(nCase,1);
ValidCase = false(nCase,1);

case_masks = false(rows,cols,nCase);
case_id = 0;

for ii = 1:numel(initial_dilation_km_list)
    init_km = initial_dilation_km_list(ii);
    init_px = round(init_km*1000/pixel_res);

    init_mask = imdilate(seed_hull_c, ...
        strel('disk',init_px,0));

    init_mask = init_mask & crop_roi_c;

    for ib = 1:numel(contraction_bias_list)
        bias = contraction_bias_list(ib);
        case_id = case_id+1;

        fprintf('[%d/%d] init=%.1f km, bias=%+.2f\n', ...
            case_id,nCase,init_km,bias);

        raw_c = activecontour( ...
            edge_feature_c,init_mask,n_iterations, ...
            'edge', ...
            'SmoothFactor',smooth_factor, ...
            'ContractionBias',bias);

        raw_c = raw_c & crop_roi_c;


        raw_c = keep_seed_connected(raw_c,seed_union_c);


        raw_c = imfill(raw_c,'holes');

        if close_radius_px>0
            raw_c = imclose(raw_c, ...
                strel('disk',close_radius_px,0));
        end

        if open_radius_px>0
            raw_c = imopen(raw_c, ...
                strel('disk',open_radius_px,0));
        end

        raw_c = keep_seed_connected(raw_c,seed_union_c);

        min_pixels = ceil( ...
            min_final_area_km2/pixel_area_km2);

        raw_c = bwareafilt(raw_c, ...
            [min_pixels inf],8);

        full_mask = false(rows,cols);
        full_mask(r1:r2,c1:c2) = raw_c;

        case_masks(:,:,case_id) = full_mask;

        InitialDilation_km(case_id) = init_km;
        ContractionBias(case_id) = bias;

        crop_boundary_c = bwperim(crop_roi_c,8);
        TouchesCropBoundary(case_id) = ...
            any(raw_c(:) & crop_boundary_c(:));

        if any(full_mask(:))
            G = calculate_geometry( ...
                full_mask,X_master,Y_master,pixel_res);

            Area_km2(case_id) = G.Area_km2;
            MajorAxisEllipse_km(case_id) = ...
                G.MajorAxisEllipse_km;
            MinorAxisEllipse_km(case_id) = ...
                G.MinorAxisEllipse_km;
            AlongAxisSpan_km(case_id) = ...
                G.AlongAxisSpan_km;
            AxisAzimuth_deg(case_id) = ...
                G.AxisAzimuth_deg;


            ValidCase(case_id) = ...
                ~TouchesCropBoundary(case_id) && ...
                all(full_mask(seed_union));
        end
    end
end

Sensitivity = table( ...
    InitialDilation_km,ContractionBias, ...
    Area_km2,MajorAxisEllipse_km, ...
    MinorAxisEllipse_km,AlongAxisSpan_km, ...
    AxisAzimuth_deg,TouchesCropBoundary,ValidCase);

fprintf('\n===== Edge-constrained active-contour sensitivity =====\n');
disp(Sensitivity);

valid_ids = find(ValidCase);

if numel(valid_ids)<3
    error(['Fewer than three active-contour parameter combinations are valid. ' ...
        'Inspect the sensitivity diagnostics and, if necessary, adjust ' ...
        'ContractionBias or the edge-control features.']);
end

%% 14. Consensus footprint
valid_stack = case_masks(:,:,valid_ids);
vote_fraction = mean(valid_stack,3);

ridge_footprint_consensus = ...
    vote_fraction>=consensus_main_fraction;

ridge_footprint_inner = ...
    vote_fraction>=consensus_inner_fraction;

ridge_footprint_outer = ...
    vote_fraction>=consensus_outer_fraction;


ridge_footprint_consensus = ...
    finalize_consensus(ridge_footprint_consensus, ...
    seed_union,close_radius_px,open_radius_px);

ridge_footprint_inner = ...
    finalize_consensus(ridge_footprint_inner, ...
    seed_union,close_radius_px,open_radius_px);

ridge_footprint_outer = ...
    finalize_consensus(ridge_footprint_outer, ...
    seed_union,close_radius_px,open_radius_px);

GeometryMain = calculate_geometry( ...
    ridge_footprint_consensus,X_master,Y_master,pixel_res);

GeometryInner = calculate_geometry( ...
    ridge_footprint_inner,X_master,Y_master,pixel_res);

GeometryOuter = calculate_geometry( ...
    ridge_footprint_outer,X_master,Y_master,pixel_res);

fprintf('\n===== Consensus rough-ice footprint =====\n');
fprintf('Valid parameter combinations: %d/%d\n',numel(valid_ids),nCase);

fprintf('\n75%% inner consensus:\n');
disp(GeometryInner);

fprintf('50%% main consensus:\n');
disp(GeometryMain);

fprintf('25%% outer consensus:\n');
disp(GeometryOuter);


%% 15. Save manuscript footprint and diagnostics
parameters = struct();
parameters.pixel_res_m = pixel_res;
parameters.small_smooth_sigma_px = small_smooth_sigma_px;
parameters.background_sigma_px = background_sigma_px;
parameters.texture_window = texture_window;
parameters.statistics_sample_stride = statistics_sample_stride;
parameters.bright_z_threshold = bright_z_threshold;
parameters.texture_z_threshold = texture_z_threshold;
parameters.core_score_percentile_list = core_score_percentile_list;
parameters.core_score_weights = core_score_weights;
parameters.min_core_area_km2 = min_core_area_km2;
parameters.core_close_radius_px = core_close_radius_px;
parameters.chosen_core_percentile = chosen_core_percentile;
parameters.chosen_core_threshold = chosen_core_threshold;
parameters.core_link_radius_m = core_link_radius_m;
parameters.crop_margin_km = crop_margin_km;
parameters.feature_weights = feature_weights;
parameters.feature_smooth_sigma_px = feature_smooth_sigma_px;
parameters.initial_dilation_km_list = initial_dilation_km_list;
parameters.contraction_bias_list = contraction_bias_list;
parameters.n_iterations = n_iterations;
parameters.smooth_factor = smooth_factor;
parameters.close_radius_px = close_radius_px;
parameters.open_radius_px = open_radius_px;
parameters.min_final_area_km2 = min_final_area_km2;
parameters.consensus_main_fraction = consensus_main_fraction;
parameters.consensus_inner_fraction = consensus_inner_fraction;
parameters.consensus_outer_fraction = consensus_outer_fraction;
parameters.valid_case_ids = valid_ids;
parameters.definition = '50% edge-constrained active-contour consensus rough-ice footprint';

final_file = fullfile(out_dir,'ridge_edge_consensus_final.mat');
save(final_file, ...
    'ridge_footprint_consensus','ridge_footprint_inner','ridge_footprint_outer', ...
    'vote_fraction','case_masks','Sensitivity', ...
    'GeometryMain','GeometryInner','GeometryOuter', ...
    'seed_union','seed_hull','edge_feature', ...
    'P_bright_rough','P_rough','core_score', ...
    'mean_anomaly_z','mean_texture_z','mean_HH','valid_count', ...
    'CandidateTable','CoreGroupTable','X_master','Y_master','file_names','parameters','-v7.3');

writetable(CandidateTable,fullfile(out_dir,'RoughIce_Core_Candidates.csv'));
writetable(CoreGroupTable,fullfile(out_dir,'RoughIce_Core_Groups.csv'));
writetable(Sensitivity,fullfile(out_dir,'RoughIce_ActiveContour_Sensitivity.csv'));
writetable(GeometryMain,fullfile(out_dir,'RoughIce_Geometry_50pct.csv'));
writetable(GeometryInner,fullfile(out_dir,'RoughIce_Geometry_75pct.csv'));
writetable(GeometryOuter,fullfile(out_dir,'RoughIce_Geometry_25pct.csv'));

%% 16. Diagnostic overview
x_km = double(X_master(1,:))/1000;
y_km = double(Y_master(:,1))/1000;

fig = figure('Color','w','Position',[30 30 1600 900]);
subplot(2,2,1);
imagesc(x_km,y_km,mean_HH); axis image; set(gca,'YDir','normal');
colormap(gca,gray); hold on;
contour(x_km,y_km,seed_union,[0.5 0.5],'y','LineWidth',1);
contour(x_km,y_km,ridge_footprint_outer,[0.5 0.5],'b--','LineWidth',1.2);
contour(x_km,y_km,ridge_footprint_consensus,[0.5 0.5],'r','LineWidth',2);
contour(x_km,y_km,ridge_footprint_inner,[0.5 0.5],'g--','LineWidth',1.2);
title('Temporal mean HH and consensus footprint');
xlabel('Easting (km)'); ylabel('Northing (km)');

subplot(2,2,2);
imagesc(x_km,y_km,edge_feature,[0 1]); axis image; set(gca,'YDir','normal');
colorbar; hold on;
contour(x_km,y_km,ridge_footprint_consensus,[0.5 0.5],'r','LineWidth',2);
title('Edge-control field'); xlabel('Easting (km)'); ylabel('Northing (km)');

subplot(2,2,3);
imagesc(x_km,y_km,vote_fraction,[0 1]); axis image; set(gca,'YDir','normal');
colorbar; hold on;
contour(x_km,y_km,ridge_footprint_consensus,[0.5 0.5],'r','LineWidth',2);
title('Fraction of valid parameter cases'); xlabel('Easting (km)'); ylabel('Northing (km)');

subplot(2,2,4);
imagesc(x_km,y_km,core_score); axis image; set(gca,'YDir','normal');
colorbar; hold on;
contour(x_km,y_km,seed_union,[0.5 0.5],'y','LineWidth',1.2);
contour(x_km,y_km,ridge_footprint_consensus,[0.5 0.5],'r','LineWidth',2);
title(sprintf('Core score and final 50%% footprint (P%.1f)',chosen_core_percentile));
xlabel('Easting (km)'); ylabel('Northing (km)');

saveas(fig,fullfile(out_dir,'RoughIce_Footprint_Overview.png'));

fig2 = figure('Color','w','Position',[30 30 1500 1100]);
for i = 1:nCase
    subplot(numel(initial_dilation_km_list),numel(contraction_bias_list),i);
    imagesc(x_km,y_km,mean_HH); axis image off; set(gca,'YDir','normal');
    colormap(gca,gray); hold on;
    if any(case_masks(:,:,i),'all')
        contour(x_km,y_km,case_masks(:,:,i),[0.5 0.5],'r','LineWidth',1.2);
    end
    title(sprintf('init %.1f km, bias %+.2f\nA=%.1f km2, valid=%d', ...
        Sensitivity.InitialDilation_km(i),Sensitivity.ContractionBias(i), ...
        Sensitivity.Area_km2(i),Sensitivity.ValidCase(i)));
end
saveas(fig2,fullfile(out_dir,'RoughIce_ActiveContour_Cases.png'));

fprintf('\n===== Rough-ice footprint complete =====\n');
fprintf('Final 50%% consensus area: %.3f km2\n',GeometryMain.Area_km2);
fprintf('Final file: %s\n',final_file);

%% Local functions
function [center,scale] = robust_center_scale(x)
    x = x(isfinite(x));

    if isempty(x)
        center = 0;
        scale = 1;
        return;
    end

    center = median(x);
    mad_value = median(abs(x-center));
    scale = 1.4826*mad_value;

    if ~isfinite(scale) || scale<1e-6
        scale = std(x,0,'omitnan');
    end

    if ~isfinite(scale) || scale<1e-6
        scale = 1;
    end
end


function out = robust_rescale01(in,roi,low_p,high_p)
    in = double(in);
    valid = roi & isfinite(in);
    values = in(valid);

    if isempty(values)
        out = zeros(size(in));
        return;
    end

    low_value = prctile(values,low_p);
    high_value = prctile(values,high_p);

    if ~isfinite(low_value) || ...
            ~isfinite(high_value) || ...
            high_value<=low_value
        low_value = min(values);
        high_value = max(values);
    end

    if high_value<=low_value
        out = zeros(size(in));
        return;
    end

    out = (in-low_value)/(high_value-low_value);
    out = max(0,min(1,out));
    out(~isfinite(out)) = 0;
end

function mask = keep_seed_connected(mask,seed)
    CC = bwconncomp(mask,8);
    keep = false(size(mask));

    for k = 1:CC.NumObjects
        pix = CC.PixelIdxList{k};

        if any(seed(pix))
            keep(pix) = true;
        end
    end

    mask = keep;
end

function mask = finalize_consensus( ...
    mask,seed,close_radius_px,open_radius_px)

    mask = mask | seed;
    mask = keep_seed_connected(mask,seed);
    mask = imfill(mask,'holes');

    if close_radius_px>0
        mask = imclose(mask, ...
            strel('disk',close_radius_px,0));
    end

    if open_radius_px>0
        mask = imopen(mask, ...
            strel('disk',open_radius_px,0));
    end

    mask = keep_seed_connected(mask,seed);
end

function Geometry = calculate_geometry(mask,X,Y,pixel_res)
    x = double(X(mask));
    y = double(Y(mask));

    x_mean = mean(x);
    y_mean = mean(y);

    x0 = x-x_mean;
    y0 = y-y_mean;

    C = cov([x0 y0]);
    [V,D] = eig(C);
    [~,principal_id] = max(diag(D));
    principal_vector = V(:,principal_id);

    along = ...
        x0*principal_vector(1)+ ...
        y0*principal_vector(2);

    across = ...
        -x0*principal_vector(2)+ ...
         y0*principal_vector(1);

    axis_azimuth = mod( ...
        atan2d(principal_vector(1), ...
               principal_vector(2)),180);

    major_axis = ...
        4*sqrt(max(var(along,1),0));

    minor_axis = ...
        4*sqrt(max(var(across,1),0));

    Geometry = table( ...
        nnz(mask)*pixel_res^2/1e6, ...
        major_axis/1000, ...
        minor_axis/1000, ...
        range(along)/1000, ...
        x_mean,y_mean,axis_azimuth, ...
        'VariableNames',{ ...
        'Area_km2','MajorAxisEllipse_km', ...
        'MinorAxisEllipse_km','AlongAxisSpan_km', ...
        'CentroidX_m','CentroidY_m', ...
        'AxisAzimuth_deg'});
end
