%% 02_train_validate_random_forest.m
% Train and validate the random-forest polynya classifier used in the
% ESSiP manuscript.
%
% Thirty manually labelled spatiotemporal patches (123,516 labelled pixels)
% are assigned intact to training and validation groups. After validation,
% the final classifier is refit using all labelled patches for scene-wide
% classification.
%
% Manuscript settings:
%   - predictors: Gaussian-smoothed HH backscatter and 5 x 5 local HH SD
%   - number of trees: 300
%   - minimum leaf size: 5
%   - final classification probability threshold: 0.60
%   - fixed random seeds are used for reproducibility
%
% MATLAB R2022b
% Required toolboxes:
%   - Image Processing Toolbox
%   - Statistics and Machine Learning Toolbox

clear; clc; close all;


%% Repository paths
project_root = fileparts(mfilename('fullpath'));


%% 1. Configuration
train_dir = fullfile(project_root,'training_patches');
model_dir = fullfile(project_root,'model');

model_out   = fullfile(model_dir,'RF_model_clean_final.mat');
report_csv  = fullfile(model_dir,'RF_validation_metrics_clean.csv');
loso_csv    = fullfile(model_dir,'RF_LOSO_per_scene_clean.csv');
manifest_csv = fullfile(model_dir,'RF_patch_manifest_clean.csv');

if ~exist(model_dir,'dir'), mkdir(model_dir); end

scene_tags = { ...
    '20260223','20260323','20251213','20251219', ...
    '20260109','20260122','20260404','20260412', ...
    '20260416','20260425'};

feature_names = {'HH_s','HH_std'};

rng_split = 1;
rng_train = 42;
nTrees = 300;
minLeaf = 5;

th_final = 0.60;
th_list = [0.30 0.50 0.60 0.70];

run_LOSO = true;
save_cv_model = false;

expected_patches = 30;
expected_samples = 123516;

%% 2. Load all labelled training patches
patchData = struct('scene',{},'patchID',{},'X',{},'y',{},'nSamples',{});
patch_counter = 0;

for f = 1:numel(scene_tags)
    tag = scene_tags{f};

    mask_file = fullfile(train_dir,['TrainPatches_' tag '_mask.mat']);
    crop_file = fullfile(train_dir,['Cropped_' tag '_patches.mat']);

    if ~isfile(mask_file), error('Missing label file: %s',mask_file); end
    if ~isfile(crop_file), error('Missing patch-image file: %s',crop_file); end

    A = load(mask_file,'TrainingPatches');
    B = load(crop_file,'patch_HH');

    TrainingPatches = A.TrainingPatches;
    patch_HH = B.patch_HH;

    if numel(TrainingPatches) ~= numel(patch_HH)
        error('%s: TrainingPatches and patch_HH contain different numbers of patches.',tag);
    end

    n_loaded_this_scene = 0;

    for k = 1:numel(TrainingPatches)
        if ~isfield(TrainingPatches(k),'label') || isempty(TrainingPatches(k).label)
            continue;
        end

        HH = double(patch_HH{k});
        label = double(TrainingPatches(k).label);

        if ~isequal(size(HH),size(label))
            error('%s patch %d: HH image and label have different sizes.',tag,k);
        end

        HH_s = imgaussfilt(HH,1,'Padding','replicate');
        HH_std = stdfilt(HH,ones(5,5));

        X_raw = [HH_s(:),HH_std(:)];
        y_raw = label(:);

        valid = isfinite(X_raw(:,1)) & isfinite(X_raw(:,2)) & ...
                (y_raw==0 | y_raw==1);

        if ~any(valid)
            warning('%s patch %d has no valid labelled samples and was skipped.',tag,k);
            continue;
        end

        patch_counter = patch_counter+1;
        n_loaded_this_scene = n_loaded_this_scene+1;

        patchData(patch_counter).scene = tag;
        patchData(patch_counter).patchID = k;
        patchData(patch_counter).X = X_raw(valid,:);
        patchData(patch_counter).y = y_raw(valid);
        patchData(patch_counter).nSamples = sum(valid);
    end

    fprintf('Loaded %s: %d valid patches\n',tag,n_loaded_this_scene);
end

if isempty(patchData), error('No labelled training patches were loaded.'); end

X_all = vertcat(patchData.X);
y_all = vertcat(patchData.y);

fprintf('\n===== All labelled samples =====\n');
fprintf('Patches: %d\n',numel(patchData));
fprintf('Samples: %d\n',size(X_all,1));
fprintf('Predictors: %d\n',size(X_all,2));
fprintf('Class 0: %d\n',sum(y_all==0));
fprintf('Class 1: %d\n',sum(y_all==1));

if numel(patchData) ~= expected_patches
    warning('Loaded %d patches; manuscript dataset contains %d.',numel(patchData),expected_patches);
end
if size(X_all,1) ~= expected_samples
    warning('Loaded %d valid samples; manuscript dataset contains %d.',size(X_all,1),expected_samples);
end

%% 3. Fixed 21/9 patch-level holdout split
rng(rng_split,'twister');

nPatch = numel(patchData);
idx = randperm(nPatch);
nVal = max(1,round(0.30*nPatch));

valPatchIdx = idx(1:nVal);
trainPatchIdx = idx(nVal+1:end);

X_train = vertcat(patchData(trainPatchIdx).X);
y_train = vertcat(patchData(trainPatchIdx).y);
X_val = vertcat(patchData(valPatchIdx).X);
y_val = vertcat(patchData(valPatchIdx).y);

fprintf('\n===== Patch-level holdout split =====\n');
fprintf('Training patches: %d; samples: %d\n',numel(trainPatchIdx),size(X_train,1));
fprintf('Validation patches: %d; samples: %d\n',numel(valPatchIdx),size(X_val,1));

%% 4. Train the 21-patch validation model
fprintf('\n===== Train 21-patch validation model =====\n');
rng(rng_train,'twister');

RF_cv = TreeBagger(nTrees,X_train,y_train, ...
    'Method','classification', ...
    'MinLeafSize',minLeaf, ...
    'OOBPrediction','On', ...
    'OOBPredictorImportance','On');

p1_val = predict_class1_probability(RF_cv,X_val);

metric_holdout = table( ...
    'Size',[numel(th_list),10], ...
    'VariableTypes',repmat({'double'},1,10), ...
    'VariableNames',{'Threshold','Accuracy','Precision','Recall','F1', ...
                     'IoU','TN','FP','FN','TP'});

for i = 1:numel(th_list)
    metric_holdout(i,:) = metric_row(y_val,p1_val,th_list(i));
end

disp(metric_holdout);

[~,best_idx] = max(metric_holdout.F1);
th_best_holdout = metric_holdout.Threshold(best_idx);

fprintf('Candidate threshold with highest holdout F1: %.2f\n',th_best_holdout);
fprintf('Fixed manuscript threshold: %.2f\n',th_final);

row_final = find(abs(metric_holdout.Threshold-th_final)<1e-12,1);
if isempty(row_final), error('The fixed threshold %.2f is absent from th_list.',th_final); end
metrics_patch_holdout = metric_holdout(row_final,:);

%% 5. Optional leave-one-scene-out (LOSO) validation
metrics_LOSO_overall = table();
metrics_LOSO_by_scene = table();
LOSO_true = [];
LOSO_prob = [];

if run_LOSO
    fprintf('\n===== Leave-one-scene-out validation (LOSO) =====\n');

    unique_scenes = unique(string({patchData.scene}),'stable');
    nScene = numel(unique_scenes);

    metrics_LOSO_by_scene = table( ...
        'Size',[nScene,11], ...
        'VariableTypes',{'string','double','double','double','double', ...
                         'double','double','double','double','double','double'}, ...
        'VariableNames',{'Scene','NSamples','Accuracy','Precision','Recall', ...
                         'F1','IoU','TN','FP','FN','TP'});

    for s = 1:nScene
        test_scene = unique_scenes(s);
        is_test_patch = string({patchData.scene}) == test_scene;

        X_fold_train = vertcat(patchData(~is_test_patch).X);
        y_fold_train = vertcat(patchData(~is_test_patch).y);
        X_fold_test = vertcat(patchData(is_test_patch).X);
        y_fold_test = vertcat(patchData(is_test_patch).y);

        fprintf('[%d/%d] Hold out scene %s: train %d, test %d\n', ...
            s,nScene,test_scene,size(X_fold_train,1),size(X_fold_test,1));

        rng(rng_train+s,'twister');
        RF_fold = TreeBagger(nTrees,X_fold_train,y_fold_train, ...
            'Method','classification', ...
            'MinLeafSize',minLeaf, ...
            'OOBPrediction','Off');

        p_fold = predict_class1_probability(RF_fold,X_fold_test);
        row = metric_row(y_fold_test,p_fold,th_final);

        metrics_LOSO_by_scene.Scene(s) = test_scene;
        metrics_LOSO_by_scene.NSamples(s) = numel(y_fold_test);
        metrics_LOSO_by_scene{s,3:end} = row{1,2:end};

        LOSO_true = [LOSO_true;y_fold_test]; %#ok<AGROW>
        LOSO_prob = [LOSO_prob;p_fold]; %#ok<AGROW>

        clear RF_fold X_fold_train y_fold_train X_fold_test y_fold_test p_fold
    end

    metrics_LOSO_overall = metric_row(LOSO_true,LOSO_prob,th_final);
    metrics_LOSO_overall.Properties.RowNames = {'LOSO_overall'};

    fprintf('\nOverall LOSO performance:\n');
    disp(metrics_LOSO_overall);
    writetable(metrics_LOSO_by_scene,loso_csv);
end

%% 6. Refit the final mapping model using all 30 patches
fprintf('\n===== Refit final mapping model using all labelled patches =====\n');
rng(rng_train,'twister');

RF_final = TreeBagger(nTrees,X_all,y_all, ...
    'Method','classification', ...
    'MinLeafSize',minLeaf, ...
    'OOBPrediction','On', ...
    'OOBPredictorImportance','On');

classNames = RF_final.ClassNames;
importance_final = RF_final.OOBPermutedPredictorDeltaError;

fprintf('RF_final training samples: %d\n',size(RF_final.X,1));
fprintf('RF_final trees: %d\n',RF_final.NumTrees);

%% 7. Save model, validation results, and provenance
PatchIndex = (1:nPatch)';
Scene = string({patchData.scene})';
PatchID = [patchData.patchID]';
NSamples = [patchData.nSamples]';
Split = repmat("train",nPatch,1);
Split(valPatchIdx) = "validation";

patch_manifest = table(PatchIndex,Scene,PatchID,NSamples,Split);
writetable(patch_manifest,manifest_csv);

validation_metrics = metric_holdout;
validation_metrics.Scheme = repmat("Patch_holdout",height(validation_metrics),1);

if run_LOSO
    loso_for_csv = metrics_LOSO_overall;
    loso_for_csv.Scheme = "LOSO_overall";
    validation_metrics = [validation_metrics;loso_for_csv];
end
writetable(validation_metrics,report_csv);

provenance = struct();
provenance.created_utc = datetime('now','TimeZone','UTC');
provenance.script = mfilename;
provenance.train_dir = train_dir;
provenance.scene_tags = scene_tags;
provenance.feature_names = feature_names;
provenance.rng_split = rng_split;
provenance.rng_train = rng_train;
provenance.nTrees = nTrees;
provenance.minLeaf = minLeaf;
provenance.th_final = th_final;
provenance.th_list = th_list;
provenance.nPatches = nPatch;
provenance.nSamples = size(X_all,1);
provenance.nClass0 = sum(y_all==0);
provenance.nClass1 = sum(y_all==1);
provenance.trainPatchIdx = trainPatchIdx;
provenance.valPatchIdx = valPatchIdx;
provenance.run_LOSO = run_LOSO;

fprintf('\nSaving final model...\n');
save(model_out, ...
    'RF_final','th_final','nTrees','minLeaf','feature_names','classNames', ...
    'importance_final','metric_holdout','metrics_patch_holdout', ...
    'metrics_LOSO_overall','metrics_LOSO_by_scene', ...
    'patch_manifest','provenance','-v7.3');

if save_cv_model
    cv_out = fullfile(model_dir,'RF_model_clean_patch_holdout_only.mat');
    save(cv_out,'RF_cv','metric_holdout','metrics_patch_holdout', ...
        'trainPatchIdx','valPatchIdx','rng_split','rng_train', ...
        'nTrees','minLeaf','th_final','feature_names','-v7.3');
end

fprintf('\n===== Completed =====\n');
fprintf('Final model: %s\n',model_out);
fprintf('Validation summary: %s\n',report_csv);
fprintf('Patch manifest: %s\n',manifest_csv);
if run_LOSO, fprintf('LOSO per-scene results: %s\n',loso_csv); end

%% Local functions
function p1 = predict_class1_probability(RF,X)
    [~,score] = predict(RF,X);
    class_names = RF.ClassNames;
    if iscell(class_names)
        class_num = str2double(class_names);
    elseif isstring(class_names)
        class_num = str2double(cellstr(class_names));
    else
        class_num = double(class_names);
    end
    col1 = find(class_num==1,1);
    if isempty(col1), error('Class 1 was not found in RF.ClassNames.'); end
    p1 = score(:,col1);
end

function row = metric_row(y_true,p1,threshold)
    y_true = double(y_true(:));
    y_pred = double(p1(:)>=threshold);

    C = confusionmat(y_true,y_pred,'Order',[0 1]);
    TN = C(1,1); FP = C(1,2);
    FN = C(2,1); TP = C(2,2);

    accuracy = (TP+TN)/max(sum(C(:)),1);
    precision = TP/max(TP+FP,1);
    recall = TP/max(TP+FN,1);
    F1 = 2*precision*recall/max(precision+recall,eps);
    IoU = TP/max(TP+FP+FN,1);

    row = table(threshold,accuracy,precision,recall,F1,IoU,TN,FP,FN,TP, ...
        'VariableNames',{'Threshold','Accuracy','Precision','Recall','F1', ...
                         'IoU','TN','FP','FN','TP'});
end
