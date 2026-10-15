%% 07_plot_figure1.m
% Purpose: Reproduce manuscript Figure 1.
% Figure 1 shows the study area, correspondence between the grounded ice ridge
% and seafloor topography, and ESSiP area evolution.
% Figure 1(a) requires M_Map, which is not redistributed with this repository.
% MATLAB: R2022b

clear; clc; close all;

%% ==================== 0. Paths and common figure settings ====================
% Public-release version for manuscript Figure 1.
% The plotting logic is unchanged; machine-specific paths are replaced by
% repository-relative paths. Figure 1(a) uses M_Map (not redistributed here).

project_root = fileparts(mfilename('fullpath'));
processed_dir = fullfile(project_root,'processed');
heat_dir = fullfile(processed_dir,'SAR_constrained_heat_budget_ERA5_cells_expertQC');

% Derived SAR background archived with the companion data record. Required
% Historical SAR background used for manuscript Figure 1(b).
% Required fields: X_master, Y_master, valid_mask, temporal_median.
sar_bg_file = fullfile(processed_dir,'ridge','iceberg_mask_for3D.mat');
ridge_file = fullfile(processed_dir,'ridge','edge_consensus_footprint', ...
    'ridge_edge_consensus_final.mat');
area_file = fullfile(heat_dir,'SAR_Polynya_Observation_Area.csv');
ibcao_file = fullfile(project_root,'IBCAO','IBCAO_v4_2_400m_ice.nc');

out_dir = fullfile(project_root,'figures_output');
if ~isfolder(out_dir), mkdir(out_dir); end

assert(isfile(sar_bg_file), 'Figure 1 SAR background file not found: %s', sar_bg_file);
assert(isfile(ridge_file), 'Rough-ice-footprint file not found: %s', ridge_file);
assert(isfile(area_file), 'Polynya-area file not found: %s', area_file);
assert(isfile(ibcao_file), 'IBCAO file not found: %s', ibcao_file);
assert(exist('m_proj','file') == 2, ...
    'M_Map is required for Figure 1(a) but was not found on the MATLAB path.');

%% ==================== 0b. Common figure settings ====================
% Practical full-width GRL/AGU layout: one cohesive multi-panel figure.
font_name = 'Arial';

% Font sizes are chosen for readability at final publication size.
panel_font_size      = 14.0;
axis_font_size       = 9.5;
label_font_size      = 10.5;
annotation_font_size = 9.0;
map_main_font_size   = 11.0;
map_minor_font_size  = 9.0;
threeD_font_size     = 9.5;
colorbar_font_size   = 9.0;


nc_file = ibcao_file;


lat_lim = [69, 77];
lon_lim = [140, 179];


box_lat = [72.5, 72.9, 72.9, 72.5, 72.5];
box_lon = [159.5, 159.5, 163.0, 163.0, 159.5];


mstruct = defaultm('stereo');
mstruct.origin = [90, 0, 0];
mstruct.geoid = wgs84Ellipsoid('meters');
mstruct.truescale = 75;
mstruct = defaultm(mstruct);

x = ncread(nc_file, 'x');
y = ncread(nc_file, 'y');



lats_boundary = linspace(lat_lim(1), lat_lim(2), 50);
lons_boundary = linspace(lon_lim(1), lon_lim(2), 50);
[Lons_b_grid, Lats_b_grid] = meshgrid(lons_boundary, lats_boundary);
[x_dense, y_dense] = mfwdtran(mstruct, Lats_b_grid, Lons_b_grid);

buffer = 50000; % 50 km buffer
idx_x = find(x >= min(x_dense(:))-buffer & x <= max(x_dense(:))+buffer);
idx_y = find(y >= min(y_dense(:))-buffer & y <= max(y_dense(:))+buffer);
% ====================================================================================


step = 3; % bathymetry subsampling stride; use 2 for higher detail
start_x = idx_x(1); start_y = idx_y(1);
x_sub = x(start_x : step : idx_x(end));
y_sub = y(start_y : step : idx_y(end));

Z_subset = ncread(nc_file, 'z', [start_x, start_y], [length(x_sub), length(y_sub)], [step, step]);
Z_grid = double(Z_subset'); 


Z_grid(Z_grid > 0) = 10; 

[X_grid, Y_grid] = meshgrid(x_sub, y_sub);
[Lat_grid, Lon_grid] = minvtran(mstruct, X_grid, Y_grid);



Lon_grid(Lon_grid < 0) = Lon_grid(Lon_grid < 0) + 360;
% ==========================================================================


fig = figure('Units','centimeters', ...
    'Position',[2 2 25 15], ...
    'Color','w', ...
    'Renderer','opengl');


ax1 = axes('Parent',fig,'Position',[0.055 0.450 0.420 0.420]);
ax2 = axes('Parent',fig,'Position',[0.535 0.450 0.420 0.50]);
ax3 = axes('Parent',fig,'Position',[0.1 0.085 0.875 0.300]);


axes(ax1);
m_proj('lambert', 'long', lon_lim, 'lat', lat_lim);


m_pcolor(Lon_grid, Lat_grid, Z_grid);
shading interp; % smooth visual transitions between grid cells
hold on;


m_contour(Lon_grid, Lat_grid, Z_grid, [0, 0], 'Color', [0.1 0.1 0.1], 'LineWidth', 0.5);


m_contour(Lon_grid, Lat_grid, Z_grid, [-50:10:-10], ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 0.2, 'ShowText', 'off');



caxis([-300, 10]); 



deep_cmap = [linspace(0.0, 0.0, 500)', linspace(0.01, 0.2, 500)', linspace(0.05, 0.6, 500)'];



shallow_cmap = [linspace(0.0, 0.85, 100)', linspace(0.2, 0.95, 100)', linspace(0.6, 1.0, 100)'];


land_cmap  = repmat([0.95, 0.9, 0.8], 20, 1); 


bathy_cmap = [deep_cmap; shallow_cmap; land_cmap];
colormap(ax1,bathy_cmap);
% ===================================================================================

m_line(box_lon, box_lat, 'Color', 'r', 'LineWidth', 2.5);


m_grid('box', 'fancy', 'tickdir', 'in', 'LineWidth', 1.2, ...
    'FontSize', axis_font_size+1.5, 'FontName', font_name,'FontWeight','bold',...
    'XaxisLocation', 'bottom', 'YaxisLocation', 'left', ...
    'linestyle', 'none');

grid_text = findobj(ax1,'Type','text');
set(grid_text,'FontWeight','bold');



hold off;
%



sea_color = 'w';



m_text(166.0, 74, 'East Siberian Sea', 'FontSize', map_main_font_size+4, 'FontWeight', 'bold', ...
    'FontAngle', 'italic', 'Color', sea_color, 'HorizontalAlignment', 'center', 'FontName', 'Arial');


%     'FontAngle', 'italic', 'Color', sea_color, 'HorizontalAlignment', 'center', 'FontName', 'Arial');


% %     'Color', [0.15 0.15 0.15], 'HorizontalAlignment', 'center', 'FontName', 'Arial','Rotation', -15);


%     'Color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial');
%     'Color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial');



m_text(152.0, 70.3, 'S  I  B  E  R  I  A', 'FontSize', map_main_font_size+2, 'FontWeight', 'bold', ...
    'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center', 'FontName', 'Arial');
% =========================================================================

p1 = ax1.Position;
cb_x = p1(1) + 0.93*p1(3);
cb_y = p1(2) + 0.35*p1(4);
cb_w = 0.025*p1(3);
cb_h = 0.43*p1(4); 



ax_cb = axes('Parent',fig,'Position', [cb_x, cb_y, cb_w, cb_h]);




cb_deep = linspace(-300, -50.01, 200)';   % deeper-water colorbar segment
cb_shallow = linspace(-50, 0, 300)';      % shallow-shelf colorbar segment
cb_data = [cb_deep; cb_shallow];          % combine into the custom colorbar raster


imagesc(ax_cb, cb_data);
set(ax_cb, 'YDir', 'normal'); % values increase upward



colormap(ax_cb,bathy_cmap);
caxis(ax_cb, [-300, 10]); 




target_vals = [-300, -100, -30, -10];


tick_positions = interp1(cb_data, 1:500, target_vals);


set(ax_cb, 'YTick', tick_positions);
set(ax_cb, 'YTickLabel', {'-300', '-100', '-30', '-10'}, ...
    'FontName', font_name, 'FontSize', colorbar_font_size, 'FontWeight', 'bold');
    


set(ax_cb, 'XTick', []); % no x ticks on the custom colorbar
ax_cb.YAxisLocation = 'right'; % colorbar labels on the right
ax_cb.TickDir = 'in';
ax_cb.LineWidth = 1.0;
ax_cb.XColor = 'k'; ax_cb.YColor = 'k'; % black frame



land_y = cb_y + cb_h + 0.015;
land_h = 0.025;
annotation('rectangle', [cb_x, land_y, cb_w, land_h], ...
    'FaceColor', [0.95, 0.9, 0.8], 'EdgeColor', 'k', 'LineWidth', 1.0);
    
annotation('textbox', [cb_x + cb_w + 0.005, land_y, 0.05, land_h], ...
    'String', 'Land', 'EdgeColor', 'none', 'FontName', 'Arial', ...
    'FontSize', colorbar_font_size, 'FontWeight', 'bold', 'Color', 'k', ...
    'VerticalAlignment', 'middle', 'Margin', 0);



title_y = land_y + land_h + 0.012;
annotation('textbox', [cb_x - 0.05, title_y, cb_w + 0.1, 0.03], ...
    'String', 'Bathymetry(m)', 'EdgeColor', 'none', 'FontName', 'Arial', ...
    'FontSize', colorbar_font_size+0.5, 'FontWeight', 'bold', 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Margin', 0);




% The consensus file contains the final persistent rough-ice footprint.
SARData = load(sar_bg_file);
RidgeData = load(ridge_file);

required_sar = {'X_master','Y_master','valid_mask','temporal_median'};
for ii = 1:numel(required_sar)
    if ~isfield(SARData,required_sar{ii})
        error('SAR background file is missing variable: %s',required_sar{ii});
    end
end

if ~isfield(RidgeData,'X_master') || ~isfield(RidgeData,'Y_master')
    error('Rough-ice file is missing X_master/Y_master.');
end



if isfield(RidgeData,'ridge_footprint_consensus')
    ridge_mask = logical(RidgeData.ridge_footprint_consensus);
elseif isfield(RidgeData,'ridge_footprint_inner')
    ridge_mask = logical(RidgeData.ridge_footprint_inner);
else
    error('Rough-ice file is missing ridge_footprint_consensus/ridge_footprint_inner.');
end

X_master = double(SARData.X_master);
Y_master = double(SARData.Y_master);
valid_mask = logical(SARData.valid_mask);
temporal_median = double(SARData.temporal_median);


assert(isequal(size(ridge_mask),size(valid_mask)), ...
    'Rough-ice mask and SAR background have different dimensions.');
assert(isequal(size(RidgeData.X_master),size(X_master)), ...
    'Rough-ice X_master and SAR background have different dimensions.');
assert(isequal(size(RidgeData.Y_master),size(Y_master)), ...
    'Rough-ice Y_master and SAR background have different dimensions.');

grid_tol_m = 1.0;
assert(max(abs(double(RidgeData.X_master(:))-X_master(:)),[],'omitnan') < grid_tol_m, ...
    'Rough-ice and SAR-background X_master grids differ.');
assert(max(abs(double(RidgeData.Y_master(:))-Y_master(:)),[],'omitnan') < grid_tol_m, ...
    'Rough-ice and SAR-background Y_master grids differ.');

proj_sar = projcrs(3413);
[Lat_sar, Lon_sar] = projinv(proj_sar, X_master, Y_master);

idx_valid = find(valid_mask);
lon_min = min(Lon_sar(idx_valid));
lon_max = max(Lon_sar(idx_valid));
lat_min = min(Lat_sar(idx_valid));
lat_max = max(Lat_sar(idx_valid));

sar_clim = [prctile(temporal_median(valid_mask), 1), prctile(temporal_median(valid_mask), 99)];
sar_norm = (temporal_median - sar_clim(1)) / (sar_clim(2) - sar_clim(1));
sar_norm(sar_norm < 0) = 0; sar_norm(sar_norm > 1) = 1;
sar_rgb = repmat(sar_norm, [1, 1, 3]); 
% =========================================================================

% =========================================================================
nc_file = ibcao_file;
proj_bathy = projcrs(3996); 



[Lon_bathy, Lat_bathy] = meshgrid(linspace(lon_min, lon_max, 500), linspace(lat_min, lat_max, 500));
[X_reg, Y_reg] = projfwd(proj_bathy, Lat_bathy, Lon_bathy);


x_raw = ncread(nc_file, 'x'); y_raw = ncread(nc_file, 'y');
x_min = min(X_reg(:)); x_max = max(X_reg(:));
y_min = min(Y_reg(:)); y_max = max(Y_reg(:));


ix = find(x_raw >= x_min-2000 & x_raw <= x_max+2000);
iy = find(y_raw >= y_min-2000 & y_raw <= y_max+2000);

Z_raw = ncread(nc_file, 'z', [min(ix), min(iy)], [length(ix), length(iy)]);
Z_raw = double(Z_raw'); 
[X_raw_mesh, Y_raw_mesh] = meshgrid(x_raw(ix), y_raw(iy));


F = scatteredInterpolant(X_raw_mesh(:), Y_raw_mesh(:), Z_raw(:), 'natural', 'none');
Z_bathy = F(X_reg, Y_reg);

[max_z, linear_idx] = max(Z_bathy(:),[],'omitnan'); 
[row_p, col_p] = ind2sub(size(Z_bathy), linear_idx);
peak_lat = Lat_bathy(row_p, col_p);
peak_lon = Lon_bathy(row_p, col_p);
min_z = min(Z_bathy(:),[],'omitnan');


B = bwboundaries(ridge_mask);
boundary_lines = cell(length(B), 1);
for k = 1:length(B)
    b = B{k};
    boundary_lines{k} = [Lon_sar(sub2ind(size(Lon_sar), b(:,1), b(:,2))), ...
                         Lat_sar(sub2ind(size(Lat_sar), b(:,1), b(:,2)))];
end
props = regionprops(ridge_mask, 'Centroid');



hover_z = 1;                     % SAR texture surface slightly above the bathymetry
base_z = floor(min_z/5)*5 - 1;    % base elevation below the bathymetry

Z_sar_plane = repmat(hover_z, size(Lat_sar));
Z_sar_plane(~valid_mask) = NaN; 

% =========================================================================

% =========================================================================
axes(ax2);
hold(ax2,'on'); 


h_sar = surf(Lon_sar, Lat_sar, Z_sar_plane, sar_rgb, ...
     'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.95);deg_sym = char(176); 
lon_ticks = ceil(lon_min*2)/2 : 0.5 : floor(lon_max*2)/2;
lat_ticks = ceil(lat_min*5)/5 : 0.2 : floor(lat_max*5)/5;


for l = lon_ticks
    plot3([l l], [lat_min lat_max], repmat(hover_z+0.2, 1, 2), ...
          '--', 'LineWidth', 1.2, 'Color', [0.82 0.82 0.82]);

    text(l, lat_max, hover_z+0.1, sprintf('%.1f%sE', l, deg_sym), ...
         'Color', 'k', 'FontSize', threeD_font_size, 'FontWeight', 'bold', ...
         'BackgroundColor', 'none', 'EdgeColor', 'none', 'Margin', 1, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

for l = lat_ticks

    plot3([lon_min lon_max], [l l], repmat(hover_z+0.2, 1, 2), ...
          '--', 'LineWidth', 1.2, 'Color', [0.82 0.82 0.82]);
    text(lon_max, l, hover_z+0.1, sprintf('%.1f%sN', l, deg_sym), ...
         'Color', 'k', 'FontSize', threeD_font_size, 'FontWeight', 'bold', ...
         'BackgroundColor', 'none', 'EdgeColor', 'none', 'Margin', 1, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end


surf(Lon_bathy, Lat_bathy, Z_bathy, 'EdgeColor', 'none', 'FaceColor', 'interp');

[rows, cols] = size(Z_bathy);
z_S = Z_bathy(1,:); z_N = Z_bathy(end,:);
z_W = Z_bathy(:,1)'; z_E = Z_bathy(:,end)';
wall_color = [0.15 0.15 0.15]; 


surf([lon_min lon_max; lon_min lon_max], [lat_min lat_min; lat_max lat_max], ...
     repmat(base_z, 2, 2), 'FaceColor', wall_color, 'EdgeColor', 'none');
surf([Lon_bathy(1,:); Lon_bathy(1,:)], [Lat_bathy(1,:); Lat_bathy(1,:)], ...
    [z_S; repmat(base_z, 1, cols)], 'FaceColor', wall_color, 'EdgeColor', 'none'); 
surf([Lon_bathy(end,:); Lon_bathy(end,:)], [Lat_bathy(end,:); Lat_bathy(end,:)], ...
    [z_N; repmat(base_z, 1, cols)], 'FaceColor', wall_color, 'EdgeColor', 'none'); 
surf([Lon_bathy(:,1)'; Lon_bathy(:,1)'], [Lat_bathy(:,1)'; Lat_bathy(:,1)'], ...
    [z_W; repmat(base_z, 1, rows)], 'FaceColor', wall_color, 'EdgeColor', 'none'); 
surf([Lon_bathy(:,end)'; Lon_bathy(:,end)'], [Lat_bathy(:,end)'; Lat_bathy(:,end)'], ...
    [z_E; repmat(base_z, 1, rows)], 'FaceColor', wall_color, 'EdgeColor', 'none'); 



for l = lon_ticks
    plot3([l l], [lat_min lat_min], [base_z, base_z + 2], 'w-', 'LineWidth', 1.5);
    

    halign = 'center';
    if l == lon_ticks(1)
        halign = 'center'; 
    end
    text_lat = lat_min-(lat_max-lat_min)*0.006;


if l == lon_ticks(1)
    text_lat = lat_min-(lat_max-lat_min)*0.1;
end

text(ax2,l,text_lat,base_z, ...
    sprintf('%.1f%sE',l,deg_sym), ...
    'Color','k', ...
    'FontSize',threeD_font_size, ...
    'FontWeight','bold', ...
    'FontName',font_name, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top');end

for l = lat_ticks
    plot3([lon_min lon_min], [l l], [base_z, base_z + 2], 'w-', 'LineWidth', 1.5);
    text(lon_min - (lon_max-lon_min)*0.015, l, base_z, sprintf('%.1f%sN', l, deg_sym), ...
         'Color', 'k', 'FontSize', threeD_font_size, 'FontWeight', 'bold', 'FontName', font_name, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end


corner_lons = [lon_min, lon_max, lon_max, lon_min, lon_min];
corner_lats = [lat_min, lat_min, lat_max, lat_max, lat_min];
for i = 1:4 
    plot3([corner_lons(i), corner_lons(i)], [corner_lats(i), corner_lats(i)], [base_z, hover_z], ...
          ':', 'LineWidth', 1.0, 'Color', [0.6 0.6 0.6]);
end


for k = 1:length(boundary_lines)
    b_coords = boundary_lines{k};
    plot3(b_coords(:,1), b_coords(:,2), repmat(hover_z + 0.3, size(b_coords, 1), 1), 'r-', 'LineWidth', 2.0);
end
plot3(peak_lon, peak_lat, max_z+0.3 , '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'MarkerSize', 5);

if ~isempty(props)
    c_col = round(props(1).Centroid(1)); c_row = round(props(1).Centroid(2));
    ice_center_lon = Lon_sar(c_row, c_col); ice_center_lat = Lat_sar(c_row, c_col);
    
    plot3(ice_center_lon, ice_center_lat, hover_z + 0.3, 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    plot3([ice_center_lon, peak_lon], [ice_center_lat, peak_lat], [hover_z, max_z + 0.8], ...
          'r--', 'LineWidth', 1.5, 'Color', [0.9 0 0]);
end

% =========================================================================

% =========================================================================
c_purple = [0.25, 0.00, 0.45]; 
c_blue   = [0.15, 0.40, 0.85]; 
c_cyan   = [0.35, 0.80, 0.95]; 
c_lightblue  = [0.60, 0.85, 0.98]; 
cmap_custom = interp1([0, 0.4, 1], [c_purple; c_blue; c_cyan], linspace(0,1,256));
colormap(ax2,cmap_custom); 

caxis([min_z, max_z]);
xlim([lon_min, lon_max]);
ylim([lat_min, lat_max]);
zlim([base_z, hover_z + 3]); 

lat_mid = mean([lat_min, lat_max]);
phys_x = (lon_max - lon_min) * cosd(lat_mid); 
phys_y = (lat_max - lat_min);
phys_z = max(phys_x, phys_y) * 0.65; 
pbaspect([phys_x, phys_y, phys_z]); 

ax = ax2;
ax.Visible = 'off'; 
ax.Clipping = 'off';

view(ax,-28, 18); %-30


delete(findall(ax2, 'Type', 'light')); % remove existing lights

camlight(ax,-45, 20); 
lighting(ax, 'gouraud');




lighting gouraud;
material([0.35, 0.80, 0.05, 5, 0]); 


h_sar.FaceLighting = 'none';    % keep SAR texture independent of scene lighting
h_sar.AmbientStrength = 1.0;
h_sar.DiffuseStrength = 0.0;
% ====================================================

text(peak_lon + (lon_max-lon_min)*0.001, peak_lat-(lon_max-lon_min)*0.01, max_z + 2.5, sprintf('Peak: %.1f m', max_z), ...
    'Color', 'k', 'FontWeight', 'bold', 'FontSize', annotation_font_size, 'FontName', font_name, ...
    'BackgroundColor', 'none', 'EdgeColor', 'none', 'Margin', 2); 

cb = colorbar(ax2,'Location','south');
cb.Title.String = 'Bathymetry (m)';
cb.Title.Color='w';
cb.Title.FontSize = colorbar_font_size+0.5;
cb.Title.FontWeight = 'bold';
cb.FontSize = colorbar_font_size;
cb.FontWeight = 'bold';
cb.Color = 'w';
cb.LineWidth = 1.2;
cb.Ruler.TickLabelGapOffset = -3;
p2 = ax2.Position;
cb.Position = [p2(1)+0.28*p2(3), p2(2)+0.12*p2(4), 0.35*p2(3), 0.025*p2(4)];

hold off;




%% (c)

%% read
SARPolynyaObservationArea = readtable(area_file);
area_final = SARPolynyaObservationArea.TotalArea_km2;
dates = SARPolynyaObservationArea.UTC;
if ~isdatetime(dates)
    dates = datetime(dates);
end


axes(ax3);
ax = ax3;
hold(ax, 'on');


max_area = max(area_final); 
y_max = 195; 
ylim(ax, [0 y_max]);


p1_start = datetime('2025-11-16'); p1_end = datetime('2026-02-05');
p2_start = datetime('2026-02-05'); p2_end = datetime('2026-03-18');
p3_start = datetime('2026-03-18'); p3_end = datetime('2026-04-30');


fill(ax, [p1_start p1_end p1_end p1_start], [0 0 y_max y_max], ...
    [0.88 0.95 0.99], 'EdgeColor', 'none', 'FaceAlpha', 0.6);


fill(ax, [p2_start p2_end p2_end p2_start], [0 0 y_max y_max], ...
    [0.92 0.92 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    

fill(ax, [p3_start p3_end p3_end p3_start], [0 0 y_max y_max], ...
    [1.0 0.94 0.94], 'EdgeColor', 'none', 'FaceAlpha', 0.6);


header_y = y_max * 0.85;
plot(ax, [datetime('2025-11-15') datetime('2026-04-30')], [header_y header_y], ...
    '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);


text_y = y_max * 0.91;
text(ax, p1_start + caldays(30), text_y, 'Stage I: Initial Freezing', ...
    'FontSize', annotation_font_size+2, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center', 'FontName', 'Arial');
text(ax, p2_start + caldays(20), text_y, 'Stage II: Suppression', ...
    'FontSize', annotation_font_size+2, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center', 'FontName', 'Arial');
text(ax, p3_start + caldays(22), text_y, 'Stage III: Episodic Expansion', ...
    'FontSize', annotation_font_size+2, 'FontWeight', 'bold', 'Color', [0.6 0.1 0.1], ...
    'HorizontalAlignment', 'center', 'FontName', 'Arial');


line_color = [0.85 0.32 0.09]; 


fill_dates = [dates(1); dates; dates(end)];
fill_areas = [0; area_final; 0];
fill(ax, fill_dates, fill_areas, line_color, 'FaceAlpha', 0.12, 'EdgeColor', 'none');


plot(ax, dates, area_final, '-o', 'Color', line_color, 'LineWidth', 2, ...
    'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', line_color);


[~, max_idx] = max(area_final);
peak_date = dates(max_idx);


plot(ax, peak_date, max_area, 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor', '#D00000', 'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
    

text(ax, peak_date + caldays(3), max_area, ...
    sprintf('%.0f km^2\n(Apr 04)', max_area), ...
    'FontSize', annotation_font_size+2, 'FontWeight', 'bold', 'Color', '#D00000', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontName',font_name);


ylabel(ax,'Polynya Area (km^2)','FontSize',label_font_size+2,'FontWeight','bold','FontName',font_name);

x_limits = [datetime('2025-11-16') datetime('2026-04-30')];
xlim(ax, x_limits);


set(ax,'FontSize',axis_font_size+2,'FontName',font_name,'TickDir','out', ...
    'Box','off','LineWidth',1.0,'XColor','k','YColor','k', ...
    'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], 'GridLineStyle', ':');


plot(ax, x_limits, [y_max y_max], '-k', 'LineWidth', 1.2, 'HandleVisibility', 'off');

plot(ax, [x_limits(2) x_limits(2)], [0 y_max], '-k', 'LineWidth', 1.2, 'HandleVisibility', 'off');


current_ticks = ax.XTick;
eng_labels = string(current_ticks, 'MMM-yyyy', 'en_US');
ax.XTickLabel = eng_labels;
xtickangle(ax,0);
ax.XAxis.FontWeight = 'bold';

hold(ax, 'off');


p1 = ax1.Position;
p2 = ax2.Position;
p3 = ax3.Position;

annotation(fig,'textbox',[p1(1) p1(2)+p1(4)+0.006 0.05 0.025], ...
    'String','(a)','EdgeColor','none','Margin',0, ...
    'FontName',font_name,'FontSize',panel_font_size,'FontWeight','bold');

annotation(fig,'textbox',[p2(1) p2(2)+p1(4)+0.012 0.05 0.025], ...
    'String','(b)','EdgeColor','none','Margin',0, ...
    'FontName',font_name,'FontSize',panel_font_size,'FontWeight','bold');

annotation(fig,'textbox',[p1(1) p3(2)+p3(4)+0.006 0.05 0.025], ...
    'String','(c)','EdgeColor','none','Margin',0, ...
    'FontName',font_name,'FontSize',panel_font_size,'FontWeight','bold');


png_file = fullfile(out_dir,'Figure_1_Combined.png');
tif_file = fullfile(out_dir,'Figure_1_Combined.tif');
pdf_file = fullfile(out_dir,'Figure_1_Combined.pdf');

drawnow;

exportgraphics(fig,png_file, ...
    'Resolution',600, ...
    'BackgroundColor','white');

exportgraphics(fig,tif_file, ...
    'Resolution',600, ...
    'BackgroundColor','white');

exportgraphics(fig,pdf_file, ...
    'ContentType','image', ...
    'Resolution',600, ...
    'BackgroundColor','white');

fprintf('\n===== Fig. 1 saved =====\n');
fprintf('PNG: %s\n',png_file);
fprintf('TIFF: %s\n',tif_file);
fprintf('PDF: %s\n',pdf_file);
