%% 08_plot_figure2.m
% Purpose: Reproduce manuscript Figure 2.
% Figure 2 shows SAR-derived sea ice motion and ERA5 10 m winds during
% successive expansion and contraction of ESSiP during 22-26 March 2026.
% The script reads the quality-controlled output of 06_sea_ice_drift.m; it does
% not recalculate sea-ice displacement. It also writes ERA5_Wind_By_SAR_Pair.csv
% for Figure 3.
% MATLAB: R2022b

clear; clc; close all;

%% ==================== 1. Paths ====================
project_root = fileparts(mfilename('fullpath'));
processed_dir = fullfile(project_root,'processed');
out_dir = fullfile(processed_dir,'ice_drift_event_20260322_20260326');
fig_dir = fullfile(project_root,'figures_output');
if ~isfolder(fig_dir), mkdir(fig_dir); end

era5_root = fullfile(project_root,'ERA5');
inst_files = dir(fullfile(era5_root,'**','data_stream-oper_stepType-instant.nc'));
if isempty(inst_files)
    error('No ERA5 instant files were found under %s.',era5_root);
end
[~,era5_order] = sort(string({inst_files.folder}));
inst_files = inst_files(era5_order);
era5_files = arrayfun(@(d) fullfile(d.folder,d.name),inst_files,'UniformOutput',false);

all_result_file = fullfile(out_dir, ...
    'SeaIceDrift_Event_AllResults.mat');

if ~isfile(all_result_file)
    error('Sea-ice-drift result file not found: %s',all_result_file);
end

S = load(all_result_file);

required_fields = { ...
    'PairSummary','vector_tables','P', ...
    'opening_bearing_true_deg'};

for i = 1:numel(required_fields)
    if ~isfield(S,required_fields{i})
        error('Result file is missing required variable: %s',required_fields{i});
    end
end

PairSummary = S.PairSummary;
opening_bearing_true_deg = S.opening_bearing_true_deg;
ridge_file = fullfile(processed_dir,'ridge','edge_consensus_footprint', ...
    'ridge_edge_consensus_final.mat');

proj3413 = projcrs(3413);


geo_xlim = [160.0 162.8];
geo_ylim = [72.40 72.9];

n_lon = 850;
n_lat = 520;

sar_clim = [-25 0];


vector_spacing_km = 2.8;




show_scale_each_panel = true;
show_normal_each_panel = false;
show_north_arrow = false;

% opening normal
normal_arrow_length_km = 10;
normal_start_back_km = 12.5;


reference_speed_kmd = 10;


era5_lon_bounds = [160.0 164.0];
era5_lat_bounds = [72.0 73.5];


show_era5_wind_each_panel = true;
era5_reference_speed_ms = 5.0;
era5_reference_length_km = 7.0;


show_stationary_ridge_points = true;
stationary_speed_threshold_kmd = 1.0;
stationary_point_spacing_km = 1.8;
stationary_marker_size = 22;


font_name = 'Arial';
tick_font_size = 22;
axis_font_size = 18.0;
panel_font_size = 19.0;
annotation_font_size = 18.0;
era5_font_size = 14.0;
drift_scale_font_size = 14.0;


axis_line_width = 0.65;
vector_line_width = 1;
vector_head_size = 0.55;
ridge_line_width = 1.65;
normal_line_width = 1.50;
scale_line_width = 1.5;
era5_line_width = 1.8;


vector_color = [0.90 0.08 0.08];
ridge_color = [0.00 0.82 0.88];
normal_color = [0.02 0.48 0.12];
north_color = [0.35 0.35 0.35];
era5_color = [0.05 0.05 0.05];
stationary_color = [0.05 0.35 0.95];

figure_position = [20 40 1500 860];


if ~isfile(ridge_file)
    error('Rough-ice result file not found: %s',ridge_file);
end

R = load(ridge_file);

if isfield(R,'ridge_footprint_inner')
    ridge_mask_master = logical(R.ridge_footprint_inner);
elseif isfield(R,'ridge_footprint_consensus')
    ridge_mask_master = logical(R.ridge_footprint_consensus);
else
    error('Rough-ice result is missing ridge_footprint_inner/ridge_footprint_consensus.');
end

if ~isfield(R,'X_master') || ~isfield(R,'Y_master')
    error('Rough-ice result is missing X_master/Y_master.');
end

ridge_xc = mean(double(R.X_master(ridge_mask_master)),'omitnan');
ridge_yc = mean(double(R.Y_master(ridge_mask_master)),'omitnan');
[ridge_latc,ridge_lonc] = projinv(proj3413,ridge_xc,ridge_yc);


Era5Pair = read_era5_pair_mean_wind( ...
    era5_files,era5_lon_bounds,era5_lat_bounds, ...
    PairSummary,opening_bearing_true_deg);

disp('ERA5 regional-mean winds for each SAR image pair:');
disp(Era5Pair(:,{ ...
    'StartUTC','EndUTC','U10_m_s','V10_m_s', ...
    'WindSpeed_m_s','WindBearing_true_deg', ...
    'OpeningNormalWind_m_s'}));


lon_vec = linspace(geo_xlim(1),geo_xlim(2),n_lon);
lat_vec = linspace(geo_ylim(1),geo_ylim(2),n_lat);
[LONq,LATq] = meshgrid(lon_vec,lat_vec);
[Xq,Yq] = projfwd(proj3413,LATq,LONq);


fig = figure( ...
    'Color','w', ...
    'Position',figure_position, ...
    'Renderer','painters');




layout_left = 0.055;
layout_right = 0.018;
layout_bottom = 0.060;
layout_top = 0.020;
layout_hgap = 0.000;   % horizontal gap between panels
layout_vgap = 0.028;   % vertical gap between panels

tile_w = (1-layout_left-layout_right-layout_hgap)/2;
tile_h = (1-layout_bottom-layout_top-layout_vgap)/2;

tile_positions = [ ...
    layout_left, ...
    layout_bottom+tile_h+layout_vgap, ...
    tile_w,tile_h; ...
    layout_left+tile_w+layout_hgap, ...
    layout_bottom+tile_h+layout_vgap, ...
    tile_w,tile_h; ...
    layout_left, ...
    layout_bottom, ...
    tile_w,tile_h; ...
    layout_left+tile_w+layout_hgap, ...
    layout_bottom, ...
    tile_w,tile_h];

panel_letters = {'a','b','c','d'};


n_pair = height(PairSummary);
NRawRidgePoints = zeros(n_pair,1);
NRawRidgeBelowThreshold = zeros(n_pair,1);
RawRidgeMedianSpeed_kmd = nan(n_pair,1);
RawRidgeP95Speed_kmd = nan(n_pair,1);
RawRidgeStationaryFraction = nan(n_pair,1);

for p = 1:height(PairSummary)
    ax = axes(fig, ...
        'Units','normalized', ...
        'Position',tile_positions(p,:));
    hold(ax,'on');


    date1_str = datestr(PairSummary.StartUTC(p),'yyyymmdd');
    date2_str = datestr(PairSummary.EndUTC(p),'yyyymmdd');

    pair_file = fullfile(out_dir, ...
        sprintf('SeaIceDrift_Vectors_%s_%s.mat', ...
        date1_str,date2_str));

    if ~isfile(pair_file)
        error('Image-pair result file not found: %s',pair_file);
    end

    D = load(pair_file, ...
        'VectorTable','HH1','xv','yv','ridge_mask');

    T = D.VectorTable;
    HH = double(D.HH1);
    xv = double(D.xv(:)');
    yv = double(D.yv(:));
    ridge_mask = logical(D.ridge_mask);

    [xv,yv,HH] = normalize_regular_grid_local(xv,yv,HH);
    [~,~,ridge_mask] = ...
        normalize_regular_grid_local( ...
        double(D.xv(:)'),double(D.yv(:)),ridge_mask);


    F_hh = griddedInterpolant( ...
        {yv,xv},HH,'linear','none');

    F_ridge = griddedInterpolant( ...
        {yv,xv},double(ridge_mask),'nearest','none');

    HH_geo = F_hh(Yq,Xq);
    ridge_geo = F_ridge(Yq,Xq)>0.5;


    imagesc(ax,lon_vec,lat_vec,HH_geo);
    set(ax,'YDir','normal');
    colormap(ax,gray(256));
    caxis(ax,sar_clim);


    contour(ax,lon_vec,lat_vec,double(ridge_geo), ...
        [0.5 0.5], ...
        'Color',ridge_color, ...
        'LineWidth',ridge_line_width);





    raw_dt_days = T.DeltaTime_h/24;

    raw_x2_m = T.XStart_m+T.DxRaw_m;
    raw_y2_m = T.YStart_m+T.DyRaw_m;
    [raw_lat2,raw_lon2] = projinv(proj3413,raw_x2_m,raw_y2_m);
    raw_distance_m = great_circle_distance_m( ...
        T.LatStart_deg,T.LonStart_deg,raw_lat2,raw_lon2);
    raw_speed_kmd = raw_distance_m./1000./raw_dt_days;

    raw_ridge_ids = find( ...
        T.InsideStableRidge & ...
        isfinite(raw_speed_kmd));

    NRawRidgePoints(p) = numel(raw_ridge_ids);

    if ~isempty(raw_ridge_ids)
        raw_ridge_speed = raw_speed_kmd(raw_ridge_ids);
        RawRidgeMedianSpeed_kmd(p) = ...
            median(raw_ridge_speed,'omitnan');
        RawRidgeP95Speed_kmd(p) = ...
            prctile(raw_ridge_speed,95);
        RawRidgeStationaryFraction(p) = ...
            mean(raw_ridge_speed<stationary_speed_threshold_kmd);
    end

    stationary_ids = raw_ridge_ids( ...
        raw_speed_kmd(raw_ridge_ids)< ...
        stationary_speed_threshold_kmd);
    NRawRidgeBelowThreshold(p) = numel(stationary_ids);

    if show_stationary_ridge_points
        if ~isempty(stationary_ids)
            selected_stationary = spatial_thin_best_ncc( ...
                T.XStart_m(stationary_ids), ...
                T.YStart_m(stationary_ids), ...
                T.NCC(stationary_ids), ...
                stationary_point_spacing_km*1000);

            stationary_ids = stationary_ids(selected_stationary);

            scatter(ax, ...
                T.LonStart_deg(stationary_ids), ...
                T.LatStart_deg(stationary_ids), ...
                stationary_marker_size,stationary_color, ...
                'o','filled', ...
                'MarkerEdgeColor','w', ...
                'LineWidth',0.35);
        end
    end


    ids = find(T.EventAnalysisVector);

    if ~isempty(ids)
        selected = spatial_thin_best_ncc( ...
            T.XStart_m(ids),T.YStart_m(ids), ...
            T.NCC(ids),vector_spacing_km*1000);

        ids = ids(selected);


        dt_days = T.DeltaTime_h(ids)/24;

        x_end_day = T.XStart_m(ids) + ...
            T.DxCorrected_m(ids)./dt_days;

        y_end_day = T.YStart_m(ids) + ...
            T.DyCorrected_m(ids)./dt_days;

        [lat_end_day,lon_end_day] = ...
            projinv(proj3413,x_end_day,y_end_day);

        dlon = wrapTo180(lon_end_day-T.LonStart_deg(ids));
        dlat = lat_end_day-T.LatStart_deg(ids);

        quiver(ax, ...
            T.LonStart_deg(ids),T.LatStart_deg(ids), ...
            dlon,dlat,0, ...
            'Color',vector_color, ...
            'LineWidth',vector_line_width, ...
            'MaxHeadSize',vector_head_size, ...
            'AutoScale','off');
    end

%     %% 6.7 opening normal
% 
%     if show_normal_here
%         draw_opening_normal_final( ...
%             ax,ridge_latc,ridge_lonc, ...
%             opening_bearing_true_deg, ...
%             normal_arrow_length_km, ...
%             normal_color,normal_line_width, ...
%             font_name,annotation_font_size);
%     end
% 

%     if show_north_arrow && p==4
%         draw_north_arrow_final( ...
%             ax,geo_xlim,geo_ylim, ...
%             north_color,font_name,annotation_font_size);
%     end


    if show_era5_wind_each_panel
        draw_era5_wind_inset( ...
            ax,geo_xlim,geo_ylim, ...
            Era5Pair.WindSpeed_m_s(p), ...
            Era5Pair.WindBearing_true_deg(p), ...
            Era5Pair.OpeningNormalWind_m_s(p), ...
            opening_bearing_true_deg, ...
            era5_reference_speed_ms, ...
            era5_reference_length_km, ...
            era5_color,normal_color,north_color,era5_line_width, ...
            font_name,era5_font_size);
    end


    show_scale_here = show_scale_each_panel || p==4;

    if show_scale_here
        draw_velocity_reference_final( ...
            ax,geo_xlim,geo_ylim, ...
            reference_speed_kmd, ...
            vector_color,scale_line_width, ...
            font_name,drift_scale_font_size);
    end


    date_text = sprintf('%s–%s Mar 2026', ...
        datestr(PairSummary.StartUTC(p),'dd'), ...
        datestr(PairSummary.EndUTC(p),'dd'));

    panel_text = sprintf( ...
        '\\bf(%s)\\rm  %s;  U_{ice,n} = %.2f km d^{-1}', ...
        panel_letters{p},date_text, ...
        PairSummary.MedianOpeningComponent_km_day(p));

    text(ax,0.004,0.008,panel_text, ...
        'Units','normalized', ...
        'Interpreter','tex', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontName',font_name, ...
        'FontSize',panel_font_size, ...
        'FontWeight','bold', ...
        'Color','k', ...
        'BackgroundColor','w', ...
        'EdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',0.45, ...
        'Margin',2);


    xlim(ax,geo_xlim);
    ylim(ax,geo_ylim);


    mean_lat = mean(geo_ylim);
    daspect(ax,[1/cosd(mean_lat) 1 1]);

    ax.FontName = font_name;
    ax.FontSize = tick_font_size;
    ax.LineWidth = axis_line_width;
    ax.FontWeight='bold'
    ax.TickDir = 'out';
    ax.TickLength = [0.012 0.012];
    ax.Layer = 'top';
    ax.Box = 'on';
    ax.XColor = [0.12 0.12 0.12];
    ax.YColor = [0.12 0.12 0.12];

    xticks(ax,160.0:1:162.5);
    yticks(ax,72.5:0.2:72.9);

ax.XTickLabel = format_geo_ticks(ax.XTick,'E');
ax.YTickLabel = format_geo_ticks(ax.YTick,'N');


    if p<=2
        ax.XTickLabel = [];
    end

    if mod(p,2)==0
        ax.YTickLabel = [];
    end

    if p>=5
        xlabel(ax,'Longitude', ...
            'FontName',font_name, ...
            'FontSize',axis_font_size,...
            'FontWeight','bold');
    end

    if mod(p,2)==3
        ylabel(ax,'Latitude', ...
            'FontName',font_name, ...
            'FontSize',axis_font_size,...
            'FontWeight','bold');
    end
end


set(findall(fig,'-property','FontName'),'FontName',font_name);


png_file = fullfile(fig_dir, ...
    'Figure_2_SeaIceDrift.png');

pdf_file = fullfile(fig_dir, ...
    'Figure_2_SeaIceDrift.pdf');

era5_table_file = fullfile(out_dir, ...
    'ERA5_Wind_By_SAR_Pair.csv');

ridge_raw_motion_file = fullfile(out_dir, ...
    'Ridge_RawMotion_By_SAR_Pair.csv');

exportgraphics(fig,png_file, ...
    'Resolution',600, ...
    'BackgroundColor','w');

exportgraphics(fig,pdf_file, ...
    'ContentType','vector', ...
    'BackgroundColor','w');

writetable(Era5Pair,era5_table_file);

RidgeRawMotionSummary = table( ...
    PairSummary.StartUTC,PairSummary.EndUTC, ...
    NRawRidgePoints,NRawRidgeBelowThreshold, ...
    RawRidgeMedianSpeed_kmd,RawRidgeP95Speed_kmd, ...
    RawRidgeStationaryFraction, ...
    'VariableNames',{ ...
    'StartUTC','EndUTC','NQualityControlledRidgePoints', ...
    'NRawSpeedBelowThreshold','RawMedianSpeed_km_day', ...
    'RawP95Speed_km_day','RawFractionBelowThreshold'});

writetable(RidgeRawMotionSummary,ridge_raw_motion_file);

fprintf('\n===== Final manuscript figure exported =====\n');
fprintf('PNG: %s\n',png_file);
fprintf('PDF: %s\n',pdf_file);
fprintf('ERA5 statistics table: %s\n',era5_table_file);
fprintf('Uncorrected rough-ice motion statistics: %s\n',ridge_raw_motion_file);

disp('Direct uncorrected motion statistics inside the rough-ice footprint:');
disp(RidgeRawMotionSummary);



function Era5Pair = read_era5_pair_mean_wind( ...
    nc_files,lon_bounds,lat_bounds,PairSummary,opening_bearing_deg)

    if ischar(nc_files) || isstring(nc_files)
        nc_files = cellstr(nc_files);
    end

    all_time = NaT(0,1);
    all_u = [];
    all_v = [];
    lat_sub_ref = [];

    for f = 1:numel(nc_files)
        nc_file = nc_files{f};
        lon_full = double(ncread(nc_file,'longitude'));
        lat_full = double(ncread(nc_file,'latitude'));
        time_raw = double(ncread(nc_file,'valid_time'));

        lon360 = mod(lon_full,360);
        bounds360 = mod(lon_bounds,360);
        if bounds360(1)>bounds360(2)
            error('ERA5 longitude bounds crossing 0 degrees are not supported here.');
        end

        lon_idx = find(lon360>=bounds360(1) & lon360<=bounds360(2));
        lat_idx = find(lat_full>=min(lat_bounds) & lat_full<=max(lat_bounds));
        if isempty(lon_idx) || isempty(lat_idx)
            error('ERA5 file does not cover the specified study region: %s',nc_file);
        end

        start_vec = [min(lon_idx),min(lat_idx),1];
        count_vec = [numel(lon_idx),numel(lat_idx),numel(time_raw)];
        u10 = double(ncread(nc_file,'u10',start_vec,count_vec));
        v10 = double(ncread(nc_file,'v10',start_vec,count_vec));
        lat_sub = lat_full(lat_idx);

        if isempty(lat_sub_ref)
            lat_sub_ref = lat_sub;
        elseif ~isequal(size(lat_sub_ref),size(lat_sub)) || ...
                any(abs(lat_sub_ref-lat_sub)>1e-10)
            error('ERA5 latitude grids differ between input files.');
        end

        time_units = '';
        try
            time_units = ncreadatt(nc_file,'valid_time','units');
        catch
        end
        time_dt = decode_era5_time(time_raw,time_units);

        all_time = [all_time; time_dt(:)]; %#ok<AGROW>
        if isempty(all_u)
            all_u = u10;
            all_v = v10;
        else
            all_u = cat(3,all_u,u10);
            all_v = cat(3,all_v,v10);
        end
    end

    [all_time,order] = sort(all_time);
    all_u = all_u(:,:,order);
    all_v = all_v(:,:,order);
    [all_time,ia] = unique(all_time,'stable');
    all_u = all_u(:,:,ia);
    all_v = all_v(:,:,ia);

    n_pair = height(PairSummary);
    u_pair = nan(n_pair,1);
    v_pair = nan(n_pair,1);
    n_time = zeros(n_pair,1);

    for p = 1:n_pair
        pair_start = strip_datetime_timezone(PairSummary.StartUTC(p));
        pair_end = strip_datetime_timezone(PairSummary.EndUTC(p));
        tid = find(all_time>pair_start & all_time<=pair_end);
        if isempty(tid)
            error('ERA5 contains no samples for SAR interval %s to %s.', ...
                char(pair_start),char(pair_end));
        end
        u_pair(p) = weighted_space_time_mean(all_u(:,:,tid),lat_sub_ref);
        v_pair(p) = weighted_space_time_mean(all_v(:,:,tid),lat_sub_ref);
        n_time(p) = numel(tid);
    end

    speed_pair = hypot(u_pair,v_pair);
    bearing_pair = mod(atan2d(u_pair,v_pair),360);
    opening_normal_pair = ...
        u_pair.*sind(opening_bearing_deg) + ...
        v_pair.*cosd(opening_bearing_deg);

    Era5Pair = table( ...
        PairSummary.StartUTC,PairSummary.EndUTC,n_time, ...
        u_pair,v_pair,speed_pair,bearing_pair,opening_normal_pair, ...
        'VariableNames',{ ...
        'StartUTC','EndUTC','ERA5SampleCount', ...
        'U10_m_s','V10_m_s','WindSpeed_m_s', ...
        'WindBearing_true_deg','OpeningNormalWind_m_s'});
end

function time_dt = decode_era5_time(time_raw,time_units)


    units_lower = lower(strtrim(string(time_units)));
    if contains(units_lower,'hours since 1900-01-01')
        time_dt = datetime(1900,1,1) + hours(time_raw);
    elseif contains(units_lower,'hours since 1970-01-01')
        time_dt = datetime(1970,1,1) + hours(time_raw);
    else
        time_dt = datetime(1970,1,1) + seconds(time_raw);
    end
    time_dt.TimeZone = '';
end

function dt = strip_datetime_timezone(dt)
    if ~isempty(dt.TimeZone)
        dt.TimeZone = 'UTC';
        dt.TimeZone = '';
    end
end

function value = weighted_space_time_mean(A,latitudes)


    weights = reshape(cosd(double(latitudes(:))),1,[],1);
    weights = repmat(weights,size(A,1),1,size(A,3));
    valid = isfinite(A);

    if ~any(valid,'all')
        value = NaN;
        return;
    end

    value = sum(A(valid).*weights(valid)) ./ sum(weights(valid));
end

function draw_era5_wind_inset( ...
    ax,x_limits,y_limits,wind_speed,bearing_deg,normal_wind, ...
    opening_bearing_deg,reference_speed_ms,reference_length_km, ...
    color_in,normal_color_in,north_color_in,line_width,font_name,font_size)





    box_left = x_limits(2)-0.85;
    box_right = x_limits(2);%-0.025
    box_bottom = y_limits(2)-0.142;
    box_top = y_limits(2);%-0.005

    patch(ax, ...
        [box_left box_right box_right box_left], ...
        [box_bottom box_bottom box_top box_top], ...
        'w','EdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',0.45);


    north_start_lon = box_left+0.15;
    north_start_lat = box_bottom+0.085;
    [north_end_lat,north_end_lon] = destination_point( ...
        north_start_lat,north_start_lon,0,2.0);

    quiver(ax,north_start_lon,north_start_lat, ...
        wrapTo180(north_end_lon-north_start_lon), ...
        north_end_lat-north_start_lat,0, ...
        'Color',north_color_in, ...
        'LineWidth',0.85, ...
        'MaxHeadSize',1, ...
        'AutoScale','off');

    text(ax,north_end_lon,north_end_lat+0.004,'N', ...
        'FontName',font_name, ...
        'FontSize',max(8,font_size-1), ...
        'FontWeight','bold', ...
        'Color',north_color_in, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');

    center_lon = (box_left+box_right)/2;
    center_lat = box_bottom + 0.11;

    full_arrow_length_km = ...
        reference_length_km*wind_speed/reference_speed_ms;

    [lat_start,lon_start] = destination_point( ...
        center_lat,center_lon,mod(bearing_deg+180,360), ...
        full_arrow_length_km/2);
    [lat_end,lon_end] = destination_point( ...
        center_lat,center_lon,bearing_deg,full_arrow_length_km/2);

    quiver(ax,lon_start,lat_start, ...
        wrapTo180(lon_end-lon_start),lat_end-lat_start,0, ...
        'Color',color_in, ...
        'LineWidth',line_width, ...
        'MaxHeadSize',0.80, ...
        'AutoScale','off');



    if isfinite(normal_wind) && abs(normal_wind)>eps
        normal_bearing = opening_bearing_deg;
        if normal_wind<0
            normal_bearing = mod(opening_bearing_deg+180,360);
        end

        normal_arrow_length_km = ...
            reference_length_km*abs(normal_wind)/reference_speed_ms;

        [normal_lat_start,normal_lon_start] = destination_point( ...
            center_lat,center_lon,mod(normal_bearing+180,360), ...
            normal_arrow_length_km/2);
        [normal_lat_end,normal_lon_end] = destination_point( ...
            center_lat,center_lon,normal_bearing, ...
            normal_arrow_length_km/2);

        quiver(ax,normal_lon_start,normal_lat_start, ...
            wrapTo180(normal_lon_end-normal_lon_start), ...
            normal_lat_end-normal_lat_start,0, ...
            'Color',normal_color_in, ...
            'LineWidth',max(1.0,0.72*line_width), ...
            'MaxHeadSize',0.88, ...
            'AutoScale','off');
    end


    scatter(ax,center_lon,center_lat,12,color_in,'filled');

    text(ax,center_lon,box_bottom+0.040, ...
        sprintf('U_{wind} : %.1f m s^{-1}',wind_speed), ...
        'Interpreter','tex', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','bold', ...
        'Color',color_in, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');

    text(ax,center_lon,box_bottom+0.005, ...
        sprintf('U_{wind,n} : %+.2f m s^{-1}',normal_wind), ...
        'Interpreter','tex', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','bold', ...
        'Color',normal_color_in, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');
end

function selected_local_ids = spatial_thin_best_ncc(x,y,ncc,spacing_m)



    x = double(x(:));
    y = double(y(:));
    ncc = double(ncc(:));

    x0 = min(x);
    y0 = min(y);

    ix = floor((x-x0)/spacing_m);
    iy = floor((y-y0)/spacing_m);

    key = ix + 1i*iy;
    [unique_key,~,group_id] = unique(key,'stable'); 

    selected_local_ids = zeros(numel(unique_key),1);

    for g = 1:numel(unique_key)
        members = find(group_id==g);
        [~,best] = max(ncc(members));
        selected_local_ids(g) = members(best);
    end

    selected_local_ids = sort(selected_local_ids);
end

function draw_opening_normal_final( ...
    ax,ridge_lat,ridge_lon,bearing_deg,length_km, ...
    color_in,line_width,font_name,font_size)




    opposite_bearing = mod(bearing_deg+180,360);

    [lat_start,lon_start] = destination_point( ...
        ridge_lat,ridge_lon,opposite_bearing,13.5);

    [lat_end,lon_end] = destination_point( ...
        lat_start,lon_start,bearing_deg,length_km);

    quiver(ax,lon_start,lat_start, ...
        wrapTo180(lon_end-lon_start),lat_end-lat_start,0, ...
        'Color',color_in, ...
        'LineWidth',line_width, ...
        'MaxHeadSize',0.78, ...
        'AutoScale','off');


    text_lon = lon_start-0.035;
    text_lat = lat_start-0.008;

    text(ax,text_lon,text_lat, ...
        sprintf('Opening normal (%.2f°)',bearing_deg), ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','bold', ...
        'Color',color_in, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'BackgroundColor','w', ...
        'EdgeColor','none', ...
        'Margin',2);
end

function draw_velocity_reference_final( ...
    ax,x_limits,y_limits,speed_kmd,color_in, ...
    line_width,font_name,font_size)




    lon_start = x_limits(2)-0.34;
    lat_start = y_limits(1)+0.05;

    [lat_end,lon_end] = destination_point( ...
        lat_start,lon_start,90,speed_kmd);

    quiver(ax,lon_start,lat_start, ...
        wrapTo180(lon_end-lon_start),lat_end-lat_start,0, ...
        'Color',color_in, ...
        'LineWidth',line_width, ...
        'MaxHeadSize',0.70, ...
        'AutoScale','off');

    text(ax,(lon_start+lon_end)/2,lat_start-0.01, ...
        sprintf('%g km d^{-1}',speed_kmd), ...
        'Interpreter','tex', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','bold', ...
        'Color',color_in, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'BackgroundColor','w', ...
        'Margin',1);
end

function draw_north_arrow_final( ...
    ax,x_limits,y_limits,color_in,font_name,font_size)

    lon_start = x_limits(1)+0.11;
    lat_start = y_limits(1)+0.045;

    [lat_end,lon_end] = destination_point( ...
        lat_start,lon_start,0,10);

    quiver(ax,lon_start,lat_start, ...
        wrapTo180(lon_end-lon_start),lat_end-lat_start,0, ...
        'Color',color_in, ...
        'LineWidth',1.3, ...
        'MaxHeadSize',0.75, ...
        'AutoScale','off');

    text(ax,lon_end,lat_end+0.006,'N', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','bold', ...
        'Color',color_in, ...
        'HorizontalAlignment','center');
end

function [lat2,lon2] = destination_point( ...
    lat1,lon1,bearing_deg,distance_km)

    earth_radius_km = 6371.0088;

    phi1 = deg2rad(lat1);
    lambda1 = deg2rad(lon1);
    theta = deg2rad(bearing_deg);
    delta = distance_km/earth_radius_km;

    phi2 = asin( ...
        sin(phi1).*cos(delta) + ...
        cos(phi1).*sin(delta).*cos(theta));

    lambda2 = lambda1 + atan2( ...
        sin(theta).*sin(delta).*cos(phi1), ...
        cos(delta)-sin(phi1).*sin(phi2));

    lat2 = rad2deg(phi2);
    lon2 = rad2deg(lambda2);
end

function distance_m = great_circle_distance_m(lat1,lon1,lat2,lon2)


    earth_radius_m = 6371008.8;

    phi1 = deg2rad(double(lat1));
    phi2 = deg2rad(double(lat2));
    dphi = phi2-phi1;
    dlambda = deg2rad(wrapTo180(double(lon2)-double(lon1)));

    a = sin(dphi/2).^2 + ...
        cos(phi1).*cos(phi2).*sin(dlambda/2).^2;
    a = min(max(a,0),1);

    distance_m = 2*earth_radius_m*atan2(sqrt(a),sqrt(1-a));
end

function [x,y,A] = normalize_regular_grid_local(x,y,A)
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
function labels = format_geo_ticks(values,hemisphere)

labels = strings(size(values));

for i = 1:numel(values)
    if abs(values(i)-round(values(i))) < 1e-8
        labels(i) = sprintf('%.0f°%s',values(i),hemisphere);
    else
        labels(i) = sprintf('%.1f°%s',values(i),hemisphere);
    end
end
end
