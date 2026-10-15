%% 03_classify_100_scenes.m
% Purpose: Apply the final random-forest classifier to the 100-scene SAR record.
% Postprocessing removes objects smaller than 0.48 km^2 and retains complete
% candidate objects that intersect the 5 km buffer around the final
% multitemporal rough-ice footprint. Retained objects are not clipped.
% Run 01_derive_rough_ice_footprint.m and 02_train_validate_random_forest.m first.
% Inputs: acquisition manifest, processed SAR scenes, final RF model, and 50% rough-ice consensus footprint.
% MATLAB: R2022b

clear; clc; close all;


%% Repository paths
project_root = fileparts(mfilename('fullpath'));



processed_dir = fullfile(project_root,'processed');

time_csv = fullfile(processed_dir, ...
    'Sentinel1_acquisition_times_20251122_20260430.csv');

model_file = ...
    fullfile(project_root,'model','RF_model_clean_final.mat');

master_mask_file = fullfile(processed_dir, ...
    'ridge','edge_consensus_footprint', ...
    'ridge_edge_consensus_final.mat');


out_dir = fullfile(project_root,'scene_eval_RF_clean');

pixel_res = 40;
min_area_km2 = 0.48;
ridge_tolerance_m = 5000;

predict_chunk_size = 250000;
make_figures = false;  % set true for per-scene diagnostic PNGs
overwrite_existing = false;

if ~exist(out_dir,'dir'), mkdir(out_dir); end


if ~isfile(time_csv), error('Acquisition manifest not found: %s',time_csv); end

Ttime = readtable(time_csv,'TextType','string');

required_cols = ["Name","MidUTC"];
for c = required_cols
    if ~ismember(c,string(Ttime.Properties.VariableNames))
        error('Acquisition manifest is missing required column: %s',c);
    end
end

source_names = string(Ttime.Name);


if isdatetime(Ttime.MidUTC)
    Ttime.MidUTC.TimeZone = 'UTC';
else
    raw_mid = string(Ttime.MidUTC);
    parsed_mid = NaT(height(Ttime),1,'TimeZone','UTC');

    for ii = 1:height(Ttime)
        parsed_ok = false;
        candidate_formats = { ...
            'yyyy-MM-dd HH:mm:ss.S', ...
            'yyyy-MM-dd HH:mm:ss', ...
            'yyyy-MM-dd''T''HH:mm:ss.S', ...
            'yyyy-MM-dd''T''HH:mm:ss'};

        for ff = 1:numel(candidate_formats)
            try
                parsed_mid(ii) = datetime(raw_mid(ii), ...
                    'InputFormat',candidate_formats{ff}, ...
                    'TimeZone','UTC');
                parsed_ok = true;
                break;
            catch
            end
        end

        if ~parsed_ok
            error('Could not parse MidUTC at row %d: %s',ii,raw_mid(ii));
        end
    end

    Ttime.MidUTC = parsed_mid;
end

if numel(source_names) ~= 100
    warning('Acquisition manifest contains %d scenes; 100 were expected.',numel(source_names));
end
if numel(unique(source_names)) ~= numel(source_names)
    error('The acquisition manifest contains duplicate filenames.');
end

date_tokens = strings(numel(source_names),1);
for i = 1:numel(source_names)
    tok = regexp(source_names(i),'\d{8}','match','once');
    if isempty(tok), error('Could not extract a date token from filename: %s',source_names(i)); end
    date_tokens(i) = string(tok);
end

if numel(unique(date_tokens)) ~= numel(date_tokens)
    error('Multiple scenes share the same date; date-based output names would conflict.');
end

fprintf('Scenes in acquisition manifest: %d\n',numel(source_names));


M = load(model_file,'RF_final','th_final','feature_names','provenance');

if ~isfield(M,'RF_final'), error('Model file does not contain RF_final.'); end

RF = M.RF_final;
th_final = M.th_final;

fprintf('Model training samples: %d\n',size(RF.X,1));
fprintf('Number of trees: %d\n',RF.NumTrees);
fprintf('Probability threshold: %.2f\n',th_final);

if size(RF.X,2) ~= 2
    error('This classifier expects two predictors, but the model contains %d.',size(RF.X,2));
end

class_names = RF.ClassNames;
if iscell(class_names)
    class_num = str2double(class_names);
elseif isstring(class_names)
    class_num = str2double(cellstr(class_names));
else
    class_num = double(class_names);
end

col_class1 = find(class_num==1,1);
if isempty(col_class1), error('Class 1 was not found in the model.'); end


R = load(master_mask_file);

if isfield(R,'ridge_footprint_consensus')

    ridge_mask = logical(R.ridge_footprint_consensus);
    ridge_mask_source = "ridge_footprint_consensus";
elseif isfield(R,'ridge_footprint_inner')

    ridge_mask = logical(R.ridge_footprint_inner);
    ridge_mask_source = "ridge_footprint_inner";
else
    error(['Rough-ice file contains neither ridge_footprint_consensus nor ' ...
           'ridge_footprint_inner.']);
end

if ~isfield(R,'X_master') || ~isfield(R,'Y_master')
    error('Rough-ice file is missing X_master or Y_master.');
end

if ~isequal(size(ridge_mask),size(R.X_master),size(R.Y_master))
    error('Rough-ice mask dimensions do not match X_master/Y_master.');
end

fprintf('Rough-ice mask source: %s\n',ridge_mask_source);
fprintf('Rough-ice mask area (projected-pixel approximation): %.2f km2\n', ...
    nnz(ridge_mask)*(pixel_res/1000)^2);

global_dist_map = bwdist(ridge_mask)*pixel_res;
global_dist_map(global_dist_map>50000) = 50000;

min_pixels = ceil(min_area_km2*1e6/(pixel_res^2));


nScene = numel(source_names);

Date = strings(nScene,1);
SourceFile = strings(nScene,1);
MidUTC = NaT(nScene,1,'TimeZone','UTC');
ValidPixels = zeros(nScene,1);
AutoArea_km2 = nan(nScene,1);
Status = strings(nScene,1);

for i = 1:nScene
    fname = source_names(i);
    source_file = fullfile(processed_dir,fname);
    date_str = date_tokens(i);

    Date(i) = date_str;
    SourceFile(i) = fname;


    MidUTC(i) = Ttime.MidUTC(i);

    out_mat = fullfile(out_dir,['Result4_' char(date_str) '.mat']);
    out_png = fullfile(out_dir,['Eval_' char(date_str) '.png']);

    fprintf('\n[%d/%d] %s\n',i,nScene,date_str);

    if ~isfile(source_file)
        Status(i) = "missing_source";
        warning('Source file not found: %s',source_file);
        continue;
    end

    if isfile(out_mat) && ~overwrite_existing

        E = load(out_mat,'map_prob','final_auto_mask','map_classes');
        if isfield(E,'map_prob')
            ValidPixels(i) = nnz(isfinite(E.map_prob));
        end
        if isfield(E,'final_auto_mask')
            existing_mask = logical(E.final_auto_mask);
        elseif isfield(E,'map_classes')
            existing_mask = logical(E.map_classes(:,:,1));
        else
            error('Existing result is missing final_auto_mask/map_classes: %s',out_mat);
        end
        AutoArea_km2(i) = nnz(existing_mask)*(pixel_res/1000)^2;
        Status(i) = "skipped_existing";
        fprintf('Existing result found; skipping: %s\n',out_mat);
        continue;
    end

    S = load(source_file,'HH_corr_dB','xq','yq');

    if ~isfield(S,'HH_corr_dB') || ~isfield(S,'xq') || ~isfield(S,'yq')
        Status(i) = "missing_variables";
        warning('Missing HH_corr_dB/xq/yq: %s',source_file);
        continue;
    end

    HH = double(S.HH_corr_dB);
    [nr,nc] = size(HH);

    if ~isequal(size(S.xq),[nr,nc]) || ~isequal(size(S.yq),[nr,nc])
        error('%s: xq/yq dimensions do not match HH.',fname);
    end

    valid_raw = isfinite(HH);

    if nnz(valid_raw)<1000
        Status(i) = "too_few_valid_pixels";
        warning('Fewer than 1000 valid pixels; scene skipped.');
        continue;
    end


    HH_clean = HH;
    fill_value = median(HH(valid_raw),'omitnan');
    if ~isfinite(fill_value), fill_value = -15; end

    HH_clean(~valid_raw) = fill_value;
    HH_clean(HH_clean < -35) = -35;
    HH_clean(HH_clean > 0) = 0;

    HH_s = imgaussfilt(HH_clean,1,'Padding','replicate');
    HH_std = stdfilt(HH_clean,ones(5,5));


    valid_support = imerode(valid_raw,strel('square',7));

    X_full = [HH_s(:),HH_std(:)];
    valid = valid_support(:) & ...
            isfinite(X_full(:,1)) & isfinite(X_full(:,2));

    ValidPixels(i) = nnz(valid);

    if ValidPixels(i)<1000
        Status(i) = "too_few_feature_pixels";
        warning('Fewer than 1000 valid feature pixels; scene skipped.');
        continue;
    end

    p_vec = nan(nr*nc,1,'single');
    p_valid = predict_probability_chunked( ...
        RF,X_full(valid,:),col_class1,predict_chunk_size);
    p_vec(valid) = single(p_valid);
    map_prob = reshape(p_vec,nr,nc);

    raw_mask = map_prob>=th_final;
    clean_mask = bwareafilt(raw_mask,[min_pixels inf],8);

    dist_map = interp2( ...
        R.X_master,R.Y_master,double(global_dist_map), ...
        S.xq,S.yq,'linear',50000);

    CC = bwconncomp(clean_mask,8);
    final_auto_mask = false(nr,nc);

    for k = 1:CC.NumObjects
        pix = CC.PixelIdxList{k};
        if min(dist_map(pix),[],'omitnan') <= ridge_tolerance_m
            final_auto_mask(pix) = true;
        end
    end


    th_list = th_final;
    map_classes = reshape(final_auto_mask,nr,nc,1);

    postprocess_params = struct();
    postprocess_params.pixel_res_m = pixel_res;
    postprocess_params.min_area_km2 = min_area_km2;
    postprocess_params.min_pixels = min_pixels;
    postprocess_params.ridge_tolerance_m = ridge_tolerance_m;
    postprocess_params.ridge_mask_file = master_mask_file;
    postprocess_params.ridge_mask_source = ridge_mask_source;
    postprocess_params.invalid_support_erosion = '7x7 square';

    acquisition_mid_utc = MidUTC(i);
    source_name = fname;
    model_source = model_file;

    save(out_mat, ...
        'map_prob','map_classes','final_auto_mask','th_list','th_final', ...
        'acquisition_mid_utc','source_name','model_source', ...
        'postprocess_params','-v7.3');

    AutoArea_km2(i) = nnz(final_auto_mask)*(pixel_res/1000)^2;
    Status(i) = "completed";

    if make_figures
        fig = figure('Color','w','Position',[50 50 1500 520], ...
            'Visible','off');

        subplot(1,3,1);
        imagesc(HH,[-25 0]); axis image off;
        colormap(gca,gray);
        title(sprintf('%s HH',date_str));

        subplot(1,3,2);
        imagesc(map_prob,[0 1]); axis image off;
        colorbar;
        title('Class-1 probability');

        subplot(1,3,3);
        imagesc(HH,[-25 0]); axis image off;
        colormap(gca,gray); hold on;
        if any(final_auto_mask(:))
            contour(final_auto_mask,[0.5 0.5],'r','LineWidth',1.2);
        end
        title(sprintf('Final automatic mask, th=%.2f',th_final));

        try
            exportgraphics(fig,out_png,'Resolution',150);
        catch
            saveas(fig,out_png);
        end
        close(fig);
    end

    fprintf('Valid pixels: %d; automatic area: %.3f km2\n', ...
        ValidPixels(i),AutoArea_km2(i));
end


RunSummary = table(Date,SourceFile,MidUTC,ValidPixels,AutoArea_km2,Status);
summary_csv = fullfile(out_dir,'RF_clean_classification_summary.csv');
writetable(RunSummary,summary_csv);

fprintf('\n===== Classification complete =====\n');
fprintf('Completed: %d\n',sum(Status=="completed"));
fprintf('Skipped existing: %d\n',sum(Status=="skipped_existing"));
fprintf('Missing/error: %d\n',sum(~ismember(Status,["completed","skipped_existing"])));
fprintf('Run summary: %s\n',summary_csv);


function p1 = predict_probability_chunked(RF,X,col1,chunk_size)
    n = size(X,1);
    p1 = nan(n,1);
    for a = 1:chunk_size:n
        b = min(a+chunk_size-1,n);
        [~,score] = predict(RF,X(a:b,:));
        p1(a:b) = score(:,col1);
    end
end
