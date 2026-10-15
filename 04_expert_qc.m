%% 04_expert_qc.m
% Purpose: Optional interactive quality-control utility for reviewing and
% editing automatic polynya masks from 03_classify_100_scenes.m.
%
% Manual QC is observer-dependent and is not expected to reproduce the
% authors' final masks exactly when repeated independently. The released
% author-reviewed FinalMask_YYYYMMDD.mat files should be used directly for
% exact reproduction of the published downstream analyses.
%
% Input: automatic Result4_YYYYMMDD.mat files from 03_classify_100_scenes.m.
% MATLAB: R2022b
%
% Controls
%   click component : toggle component
%   A : accept all automatic candidates     C : clear selection
%   P : add polygon                         X : remove polygon
%   U : undo                                T : edit note
%   Enter/S : save and advance              Z : save confirmed zero area
%   B : previous scene                      Q/Esc : quit

clear; clc; close all;


%% Repository paths
project_root = fileparts(mfilename('fullpath'));


%% 1. Configuration
cfg = struct();

cfg.result_dir = ...
    fullfile(project_root,'scene_eval_RF_clean');

cfg.processed_dir = ...
    fullfile(project_root,'processed');

cfg.acquisition_csv = fullfile( ...
    cfg.processed_dir, ...
    'Sentinel1_acquisition_times_20251122_20260430.csv');

% The heat-budget script reads expert-confirmed masks from this subdirectory.
cfg.manual_dir = fullfile(cfg.result_dir,'manual_qc');
cfg.qc_log_csv = fullfile(cfg.manual_dir,'ManualMask_QC_Log.csv');

cfg.rf_threshold = 0.60;
cfg.display_limits_dB = [-25 -8];
cfg.resume_at_first_unreviewed = true;
cfg.start_date_token = "";  % for example, "20260322"; empty = automatic
cfg.show_component_labels = false;

if ~isfolder(cfg.result_dir)
    error('RF result directory does not exist: %s',cfg.result_dir);
end

if ~isfolder(cfg.processed_dir)
    error('Processed SAR directory does not exist: %s',cfg.processed_dir);
end

if ~isfolder(cfg.manual_dir)
    mkdir(cfg.manual_dir);
end

%% 2. Build a date-matched catalog
result_files = dir(fullfile(cfg.result_dir,'Result4_*.mat'));

if isempty(result_files)
    error('No Result4_YYYYMMDD.mat files were found in %s',cfg.result_dir);
end

date_tokens = strings(numel(result_files),1);
result_paths = strings(numel(result_files),1);

for i = 1:numel(result_files)
    date_tokens(i) = extract_date_token(result_files(i).name);
    result_paths(i) = string(fullfile( ...
        result_files(i).folder,result_files(i).name));
end

valid = strlength(date_tokens)==8;
date_tokens = date_tokens(valid);
result_paths = result_paths(valid);

[date_tokens,order] = sort(date_tokens);
result_paths = result_paths(order);

[date_tokens,ia] = unique(date_tokens,'stable');
result_paths = result_paths(ia);

Catalog = table(date_tokens,result_paths, ...
    'VariableNames',{'DateToken','ResultFile'});

Acq = read_acquisition_catalog(cfg.acquisition_csv,cfg.processed_dir);
QC = read_qc_log(cfg.qc_log_csv);

%% 3. Select the starting scene
if strlength(cfg.start_date_token)>0
    current_id = find(Catalog.DateToken==cfg.start_date_token,1);
    if isempty(current_id)
        error('Requested start date is absent: %s',char(cfg.start_date_token));
    end
elseif cfg.resume_at_first_unreviewed
    reviewed = false(height(Catalog),1);
    if ~isempty(QC)
        reviewed = ismember(Catalog.DateToken,QC.DateToken);
    end
    for i = 1:height(Catalog)
        existing_mask = fullfile(cfg.manual_dir, ...
            "FinalMask_"+Catalog.DateToken(i)+".mat");
        reviewed(i) = reviewed(i) || isfile(existing_mask);
    end
    current_id = find(~reviewed,1);
    if isempty(current_id)
        current_id = 1;
        fprintf('All %d scenes already have QC records; starting at the first scene.\n', ...
            height(Catalog));
    end
else
    current_id = 1;
end

%% 4. Interactive review loop
while current_id>=1 && current_id<=height(Catalog)
    date_token = Catalog.DateToken(current_id);
    result_file = Catalog.ResultFile(current_id);

    fprintf('\n[%d/%d] Loading %s\n', ...
        current_id,height(Catalog),char(date_token));

    Scene = load_review_scene( ...
        result_file,date_token,Acq,cfg);

    existing_file = fullfile( ...
        cfg.manual_dir,"FinalMask_"+date_token+".mat");

    [initial_mask,initial_note,was_reviewed] = ...
        load_existing_review(existing_file,Scene.candidate_mask,QC,date_token);

    action = review_one_scene( ...
        Scene,initial_mask,initial_note,was_reviewed, ...
        current_id,height(Catalog),existing_file,cfg);

    if action.Command=="save_next" || action.Command=="zero_next"
        QC = update_qc_log(QC,action.Record);
        writetable(QC,cfg.qc_log_csv);
        current_id = current_id+1;

    elseif action.Command=="previous"
        current_id = max(1,current_id-1);

    elseif action.Command=="quit"
        break;

    else
        error('Unexpected review action: %s',char(action.Command));
    end
end

%% 5. Completion summary
QC = read_qc_log(cfg.qc_log_csv);
n_reviewed = numel(intersect(Catalog.DateToken,QC.DateToken));

fprintf('\nQC progress: %d of %d scenes reviewed.\n', ...
    n_reviewed,height(Catalog));
fprintf('Manual masks: %s\n',cfg.manual_dir);
fprintf('QC log: %s\n',cfg.qc_log_csv);

if n_reviewed<height(Catalog)
    fprintf(['Review is incomplete. Do not run the final heat-budget ' ...
        'integration with automatic fallback enabled.\n']);
else
    fprintf('All RF scenes have an expert-confirmed FinalMask file.\n');
end

%% ========================================================================
% Local functions
% ========================================================================

function token = extract_date_token(text_value)
    hit = regexp(char(text_value),'\d{8}','match','once');
    if isempty(hit)
        token = "";
    else
        token = string(hit);
    end
end

function Acq = read_acquisition_catalog(csv_file,processed_dir)
    Acq = table();

    if ~isfile(csv_file)
        warning('Acquisition CSV not found; metadata and filename fallback will be used.');
        return;
    end

    T = readtable(csv_file,'TextType','string');
    vars = string(T.Properties.VariableNames);

    file_candidates = ["SourceFile","FilePath","Name","FileName","Filename"];
    time_candidates = ["MidUTC","UTC","AcquisitionMidUTC", ...
        "AcquisitionUTC","DateTime"];

    file_var = first_present(vars,file_candidates);
    time_var = first_present(vars,time_candidates);

    if strlength(file_var)==0
        warning('No recognized filename column in acquisition CSV.');
        return;
    end

    source_file = string(T.(file_var));
    for i = 1:numel(source_file)
        if ~isfile(source_file(i))
            source_file(i) = string(fullfile(processed_dir,source_file(i)));
        end
    end

    mid_utc = NaT(height(T),1,'TimeZone','UTC');
    if strlength(time_var)>0
        mid_utc = ensure_utc_datetime(T.(time_var));
    end

    date_token = strings(height(T),1);
    for i = 1:height(T)
        if ~isnat(mid_utc(i))
            date_token(i) = string(datestr(mid_utc(i),'yyyymmdd'));
        else
            date_token(i) = extract_date_token(source_file(i));
        end
    end

    Acq = table(date_token,source_file,mid_utc, ...
        'VariableNames',{'DateToken','SourceFile','MidUTC'});
end

function found = first_present(vars,candidates)
    found = "";
    for i = 1:numel(candidates)
        id = find(strcmpi(vars,candidates(i)),1);
        if ~isempty(id)
            found = vars(id);
            return;
        end
    end
end

function QC = read_qc_log(csv_file)
    if ~isfile(csv_file)
        QC = table();
        return;
    end

    QC = readtable(csv_file,'TextType','string');
    if ismember('DateToken',QC.Properties.VariableNames)
        QC.DateToken = string(QC.DateToken);
    end
    if ismember('ReviewUTC',QC.Properties.VariableNames) && ...
            ~isdatetime(QC.ReviewUTC)
        QC.ReviewUTC = ensure_utc_datetime(QC.ReviewUTC);
    end
    if ismember('Notes',QC.Properties.VariableNames)
        QC.Notes = string(QC.Notes);
        QC.Notes(ismissing(QC.Notes)) = "";
    end
end

function Scene = load_review_scene(result_file,date_token,Acq,cfg)
    info = whos('-file',char(result_file));
    names = string({info.name});

    result_source_name = "";
    result_mid_utc = NaT(1,1,'TimeZone','UTC');

    if ismember("source_name",names)
        S = load(result_file,'source_name');
        result_source_name = string(S.source_name);
        result_source_name = result_source_name(1);
    end

    if ismember("acquisition_mid_utc",names)
        S = load(result_file,'acquisition_mid_utc');
        result_mid_utc = ensure_utc_datetime(S.acquisition_mid_utc);
        result_mid_utc = result_mid_utc(1);
    end

    [candidate_mask,candidate_source,candidate_threshold] = ...
        load_candidate_mask(result_file,names,cfg.rf_threshold);

    [source_file,mid_utc] = match_source_scene( ...
        date_token,result_source_name,result_mid_utc,Acq,cfg.processed_dir);

    if ~isfile(source_file)
        error('Matched SAR source file does not exist: %s',char(source_file));
    end

    source_info = whos('-file',char(source_file));
    source_names = string({source_info.name});

    if ~ismember("HH_corr_dB",source_names)
        error('SAR source lacks HH_corr_dB: %s',char(source_file));
    end

    if all(ismember(["xv","yv"],source_names))
        S = load(source_file,'HH_corr_dB','xv','yv');
        grid_xv = double(S.xv(:)');
        grid_yv = double(S.yv(:));
    elseif all(ismember(["xq","yq"],source_names))
        S = load(source_file,'HH_corr_dB','xq','yq');
        grid_xv = double(S.xq(1,:));
        grid_yv = double(S.yq(:,1));
    else
        error('SAR source lacks xv/yv or xq/yq grid coordinates: %s',char(source_file));
    end

    HH = double(S.HH_corr_dB);
    candidate_mask = align_mask(candidate_mask,size(HH));

    Scene = struct();
    Scene.date_token = string(date_token);
    Scene.result_file = string(result_file);
    Scene.source_file = string(source_file);
    Scene.mid_utc = mid_utc;
    Scene.HH = HH;
    Scene.xv = grid_xv;
    Scene.yv = grid_yv;
    Scene.candidate_mask = logical(candidate_mask);
    Scene.candidate_source = string(candidate_source);
    Scene.candidate_threshold = candidate_threshold;
end

function [mask,source,threshold] = ...
    load_candidate_mask(result_file,names,rf_threshold)

    preferred = ["final_auto_mask","final_mask","polynya_mask"];

    for i = 1:numel(preferred)
        if ismember(preferred(i),names)
            S = load(result_file,char(preferred(i)));
            mask = logical(S.(preferred(i)));
            source = preferred(i);
            threshold = rf_threshold;
            return;
        end
    end

    if ismember("map_classes",names)
        if ismember("th_list",names)
            S = load(result_file,'map_classes','th_list');
            th_list = double(S.th_list(:));
        else
            S = load(result_file,'map_classes');
            th_list = rf_threshold;
        end

        C = move_threshold_dimension_last(S.map_classes,numel(th_list));
        [~,ith] = min(abs(th_list-rf_threshold));
        mask = logical(C(:,:,ith));
        source = "map_classes";
        threshold = th_list(ith);
        return;
    end

    error('No RF candidate mask was found in %s',char(result_file));
end

function C = move_threshold_dimension_last(C,n_threshold)
    if ndims(C)==2
        C = reshape(C,size(C,1),size(C,2),1);
        return;
    end

    s = size(C);
    threshold_dim = find(s==n_threshold,1,'last');
    if isempty(threshold_dim)
        error('Cannot identify the threshold dimension of map_classes.');
    end

    remaining = setdiff(1:ndims(C),threshold_dim,'stable');
    C = permute(C,[remaining threshold_dim]);

    if ndims(C)>3
        error('map_classes has more dimensions than expected.');
    end
end

function [source_file,mid_utc] = match_source_scene( ...
    date_token,result_source_name,result_mid_utc,Acq,processed_dir)

    source_file = "";
    mid_utc = result_mid_utc;

    if strlength(result_source_name)>0
        source_file = result_source_name;
        if ~isfile(source_file)
            source_file = string(fullfile(processed_dir,result_source_name));
        end
        if isfile(source_file)
            return;
        end
    end

    if ~isempty(Acq)
        ids = find(Acq.DateToken==date_token);
        if ~isempty(ids)
            if numel(ids)>1 && ~isnat(result_mid_utc)
                [~,jj] = min(abs(Acq.MidUTC(ids)-result_mid_utc));
                id = ids(jj);
            else
                id = ids(1);
            end
            source_file = Acq.SourceFile(id);
            if isnat(mid_utc)
                mid_utc = Acq.MidUTC(id);
            end
            if isfile(source_file)
                return;
            end
        end
    end

    candidates = dir(fullfile(processed_dir,['*' char(date_token) '*.mat']));
    if isempty(candidates)
        error(['Cannot match a SAR source scene for %s. Add source_name to ' ...
            'Result4 or correct cfg.acquisition_csv.'],char(date_token));
    end

    source_file = string(fullfile(candidates(1).folder,candidates(1).name));
end

function mask = align_mask(mask,target_size)
    mask = logical(mask);
    if isequal(size(mask),target_size)
        return;
    end

    if isequal(fliplr(size(mask)),target_size)
        mask = mask.';
        return;
    end

    error('Mask size [%s] does not match SAR size [%s].', ...
        num2str(size(mask)),num2str(target_size));
end

function [mask,note,was_reviewed] = ...
    load_existing_review(output_file,candidate_mask,QC,date_token)

    mask = false(size(candidate_mask));
    note = "";
    was_reviewed = false;

    if isfile(output_file)
        info = whos('-file',char(output_file));
        names = string({info.name});
        if ismember("final_mask",names)
            S = load(char(output_file),'final_mask');
            mask = align_mask(logical(S.final_mask),size(candidate_mask));
            was_reviewed = true;
        end
        if ismember("qc_notes",names)
            S = load(char(output_file),'qc_notes');
            note = string(S.qc_notes);
        end
    end

    if ~isempty(QC) && ismember('DateToken',QC.Properties.VariableNames)
        id = find(QC.DateToken==date_token,1,'last');
        if ~isempty(id)
            was_reviewed = true;
            if ismember('Notes',QC.Properties.VariableNames)
                note = string(QC.Notes(id));
            end
        end
    end

    % Empty cells in a CSV can be imported as MATLAB <missing> strings.
    % Normalize them to an ordinary empty string before GUI display.
    if isempty(note) || any(ismissing(note))
        note = "";
    else
        note = note(1);
    end
end

function action = review_one_scene( ...
    Scene,initial_mask,initial_note,was_reviewed, ...
    scene_id,n_scene,output_file,cfg)

    selected = logical(initial_mask);
    note = string(initial_note);
    if isempty(note) || any(ismissing(note))
        note = "";
    else
        note = note(1);
    end
    undo_stack = {};
    dirty = false;

    CC = bwconncomp(Scene.candidate_mask,8);
    L = labelmatrix(CC);

    fig = figure( ...
        'Name','RF polynya-mask expert QC', ...
        'Color','w', ...
        'NumberTitle','off', ...
        'Units','normalized', ...
        'Position',[0.04 0.06 0.92 0.86]);

    ax = axes(fig,'Position',[0.04 0.08 0.78 0.88]);

    while isvalid(fig)
        draw_review_display(ax,Scene,selected,L,CC, ...
            note,was_reviewed,scene_id,n_scene,cfg);

        figure(fig);
        [x,y,button] = ginput(1);

        if isempty(button)
            button = 13;
        end

        if button==1 || button==3
            col = round(x);
            row = round(y);
            if row>=1 && row<=size(L,1) && col>=1 && col<=size(L,2)
                component_id = L(row,col);
                if component_id>0
                    undo_stack{end+1} = selected; %#ok<AGROW>
                    pixels = CC.PixelIdxList{component_id};
                    if any(selected(pixels))
                        selected(pixels) = false;
                    else
                        selected(pixels) = true;
                    end
                    dirty = true;
                end
            end

        elseif button==double('a') || button==double('A')
            undo_stack{end+1} = selected; %#ok<AGROW>
            selected = Scene.candidate_mask;
            dirty = true;

        elseif button==double('c') || button==double('C')
            undo_stack{end+1} = selected; %#ok<AGROW>
            selected(:) = false;
            dirty = true;

        elseif button==double('p') || button==double('P')
            roi_mask = acquire_polygon_mask(ax,size(selected),[0.1 0.9 0.3]);
            if ~isempty(roi_mask)
                undo_stack{end+1} = selected; %#ok<AGROW>
                selected = selected | roi_mask;
                dirty = true;
            end

        elseif button==double('x') || button==double('X')
            roi_mask = acquire_polygon_mask(ax,size(selected),[1.0 0.3 0.2]);
            if ~isempty(roi_mask)
                undo_stack{end+1} = selected; %#ok<AGROW>
                selected(roi_mask) = false;
                dirty = true;
            end

        elseif button==double('u') || button==double('U')
            if ~isempty(undo_stack)
                selected = undo_stack{end};
                undo_stack(end) = [];
                dirty = true;
            end

        elseif button==double('t') || button==double('T')
            answer = inputdlg('QC note for this scene:', ...
                'Scene note',[3 70],{char(note)});
            if ~isempty(answer)
                note = string(answer{1});
                dirty = true;
            end

        elseif button==13 || button==double('s') || button==double('S')
            if ~any(selected(:))
                answer = questdlg( ...
                    'The selected mask is empty. Save it as confirmed zero area?', ...
                    'Confirm zero mask','Yes','No','No');
                if ~strcmp(answer,'Yes')
                    continue;
                end
            end

            Record = save_review(Scene,selected,note,output_file);
            close(fig);
            action = struct('Command',"save_next",'Record',Record);
            return;

        elseif button==double('z') || button==double('Z')
            answer = questdlg( ...
                'Confirm that the target polynya area is zero in this scene?', ...
                'Confirm zero area','Yes','No','No');
            if strcmp(answer,'Yes')
                selected(:) = false;
                Record = save_review(Scene,selected,note,output_file);
                close(fig);
                action = struct('Command',"zero_next",'Record',Record);
                return;
            end

        elseif button==double('b') || button==double('B')
            if dirty && ~confirm_discard_changes()
                continue;
            end
            close(fig);
            action = struct('Command',"previous",'Record',table());
            return;

        elseif button==27 || button==double('q') || button==double('Q')
            if dirty && ~confirm_discard_changes()
                continue;
            end
            close(fig);
            action = struct('Command',"quit",'Record',table());
            return;
        end
    end

    action = struct('Command',"quit",'Record',table());
end

function draw_review_display(ax,Scene,selected,L,CC, ...
    note,was_reviewed,scene_id,n_scene,cfg)

    cla(ax);
    imagesc(ax,Scene.HH,cfg.display_limits_dB);
    axis(ax,'image');
    axis(ax,'tight');
    set(ax,'YDir','reverse');
    colormap(ax,gray(256));
    hold(ax,'on');

    candidate_rgb = cat(3, ...
        ones(size(Scene.HH)), ...
        0.15*ones(size(Scene.HH)), ...
        0.10*ones(size(Scene.HH)));

    h_candidate = image(ax,candidate_rgb);
    set(h_candidate,'AlphaData',0.22*double(Scene.candidate_mask));

    selected_rgb = cat(3, ...
        0.10*ones(size(Scene.HH)), ...
        0.95*ones(size(Scene.HH)), ...
        0.25*ones(size(Scene.HH)));

    h_selected = image(ax,selected_rgb);
    set(h_selected,'AlphaData',0.48*double(selected));

    if any(selected(:))
        contour(ax,selected,[0.5 0.5],'Color',[0.0 1.0 0.9], ...
            'LineWidth',1.4);
    end

    if cfg.show_component_labels
        stats = regionprops(CC,'Centroid');
        for k = 1:numel(stats)
            text(ax,stats(k).Centroid(1),stats(k).Centroid(2), ...
                num2str(k),'Color','w','FontSize',7, ...
                'HorizontalAlignment','center', ...
                'BackgroundColor',[0 0 0]);
        end
    end

    review_text = "unreviewed";
    if was_reviewed
        review_text = "previously reviewed";
    end

    title(ax,{ ...
        sprintf('[%d/%d] %s | %s | candidate: %s | threshold %.2f', ...
        scene_id,n_scene,char(Scene.date_token),char(review_text), ...
        char(Scene.candidate_source),Scene.candidate_threshold), ...
        sprintf(['Candidate pixels: %d | Selected pixels: %d | Note: %s\n' ...
        'Click component: toggle | A: all | C: clear | P: add polygon | ' ...
        'X: remove polygon | U: undo | T: note | Enter: save | Z: zero | ' ...
        'B: back | Q: quit'], ...
        nnz(Scene.candidate_mask),nnz(selected),char(note))}, ...
        'Interpreter','none','FontSize',11);

    xlabel(ax,'Image column');
    ylabel(ax,'Image row');
end

function roi_mask = acquire_polygon_mask(ax,target_size,color)
    roi_mask = [];
    try
        roi = drawpolygon(ax,'Color',color,'LineWidth',1.5);
        if isempty(roi) || ~isvalid(roi)
            return;
        end
        pos = roi.Position;
        if size(pos,1)>=3
            roi_mask = poly2mask( ...
                pos(:,1),pos(:,2),target_size(1),target_size(2));
        end
        delete(roi);
    catch ME
        warning('Polygon edit cancelled or failed: %s',ME.message);
        roi_mask = [];
    end
end

function tf = confirm_discard_changes()
    answer = questdlg('Discard unsaved changes to the current scene?', ...
        'Unsaved changes','Discard','Continue editing','Continue editing');
    tf = strcmp(answer,'Discard');
end

function Record = save_review(Scene,final_mask,qc_notes,output_file)
    final_mask = logical(final_mask);
    candidate_mask = logical(Scene.candidate_mask);
    date_token = Scene.date_token;
    [~,source_base,source_ext] = fileparts(char(Scene.source_file));
    source_name = string([source_base source_ext]);
    [~,result_base,result_ext] = fileparts(char(Scene.result_file));
    rf_result_file = string([result_base result_ext]);
    acquisition_mid_utc = Scene.mid_utc;
    candidate_source = Scene.candidate_source;
    candidate_threshold = Scene.candidate_threshold;
    xv = double(Scene.xv(:)');
    yv = double(Scene.yv(:));
    qc_timestamp_utc = datetime('now','TimeZone','UTC');

    if ~any(final_mask(:))
        qc_status = "zero";
    elseif isequal(final_mask,candidate_mask)
        qc_status = "accepted_auto";
    else
        qc_status = "modified";
    end

    save(char(output_file), ...
        'final_mask','candidate_mask','date_token', ...
        'source_name','acquisition_mid_utc','rf_result_file', ...
        'candidate_source','candidate_threshold','xv','yv', ...
        'qc_status','qc_timestamp_utc','qc_notes','-v7.3');

    Record = table( ...
        date_token,qc_status, ...
        nnz(candidate_mask),nnz(final_mask), ...
        candidate_source,candidate_threshold, ...
        qc_timestamp_utc,qc_notes, ...
        string(['FinalMask_' char(date_token) '.mat']), ...
        rf_result_file,source_name, ...
        'VariableNames',{ ...
        'DateToken','Status','CandidatePixelCount','FinalPixelCount', ...
        'CandidateSource','CandidateThreshold','ReviewUTC','Notes', ...
        'OutputFile','ResultFile','SourceFile'});

    fprintf('Saved %s: %s, %d selected pixels.\n', ...
        char(date_token),char(qc_status),nnz(final_mask));
end

function QC = update_qc_log(QC,Record)
    if isempty(QC)
        QC = Record;
        return;
    end

    QC.DateToken = string(QC.DateToken);
    date_token = string(Record.DateToken);
    id = find(QC.DateToken==date_token,1);

    if isempty(id)
        QC = [QC; Record];
    else
        common = intersect( ...
            QC.Properties.VariableNames, ...
            Record.Properties.VariableNames,'stable');
        QC(id,common) = Record(:,common);
    end

    QC = sortrows(QC,'DateToken');
end

function dt = ensure_utc_datetime(raw)
    if isdatetime(raw)
        dt = raw;
        dt.TimeZone = 'UTC';
        return;
    end

    raw = string(raw);
    dt = NaT(numel(raw),1,'TimeZone','UTC');

    formats = [ ...
        "yyyy-MM-dd'T'HH:mm:ss", ...
        "yyyy-MM-dd HH:mm:ss", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSS", ...
        "yyyyMMdd'T'HHmmss"];

    for i = 1:numel(raw)
        for j = 1:numel(formats)
            try
                value = datetime(raw(i),'InputFormat',formats(j), ...
                    'TimeZone','UTC');
                if ~isnat(value)
                    dt(i) = value;
                    break;
                end
            catch
            end
        end

        if isnat(dt(i))
            try
                value = datetime(raw(i),'TimeZone','UTC');
                dt(i) = value;
            catch
            end
        end
    end
end
