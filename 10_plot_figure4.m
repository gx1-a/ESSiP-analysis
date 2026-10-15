%% 10_plot_figure4.m
% Purpose: Reproduce manuscript Figure 4.
% Figure 4 shows potential open-water heat loss, ice production, and salt
% rejection from November 2025 to April 2026. Panel (a) shows potential
% open-water net heat flux. Panel (b) shows cumulative potential ice production
% for the reference and conservative temporal-area scenarios; the right axis
% gives the corresponding potential salt rejection.
% MATLAB: R2022b

clear; close all; clc;

%% 1. Settings

project_root = fileparts(mfilename('fullpath'));
cfg.data_dir = fullfile(project_root,'processed', ...
    'SAR_constrained_heat_budget_ERA5_cells_expertQC');
cfg.out_dir = fullfile(project_root,'figures_output');

cfg.hourly_file = fullfile( ...
    cfg.data_dir,'SAR_Constrained_HeatBudget_Hourly.csv');
cfg.mat_file = fullfile( ...
    cfg.data_dir,'SAR_ERA5_Cellwise_HeatBudget.mat');

cfg.main_png = fullfile( ...
    cfg.out_dir,'Figure_4_HeatFlux_IceSalt.png');
cfg.main_pdf = fullfile( ...
    cfg.out_dir,'Figure_4_HeatFlux_IceSalt.pdf');

if ~isfolder(cfg.out_dir)
    mkdir(cfg.out_dir);
end

cfg.rho_ice_kg_m3 = 920;
cfg.water_salinity_g_kg = 30;
cfg.salt_rejection_fraction = 0.69;

cfg.main_width_cm = 20.0;
cfg.main_height_cm = 12.8;
cfg.font_name = 'Arial';
cfg.axis_font_size = 11.5;
cfg.label_font_size = 11.0;
cfg.legend_font_size = 10.0;
cfg.line_width = 1.65;
cfg.sensitivity_line_width = 1.25;
cfg.export_dpi = 600;
cfg.band_alpha = 0.20;
cfg.opening_shade_alpha = 0.09;

clr.red = [205 70 55]/255;
clr.blue = [54 116 174]/255;
clr.orange = [225 105 35]/255;
clr.grey = [145 151 157]/255;
clr.dark = [35 35 35]/255;
clr.grid = [0.86 0.87 0.88];
clr.opening = [87 160 205]/255;

%% 2. Read and validate data

assert(isfile(cfg.hourly_file), ...
    'Hourly heat-budget table not found: %s',cfg.hourly_file);
assert(isfile(cfg.mat_file), ...
    'Cellwise heat-budget MAT file not found: %s',cfg.mat_file);

H = readtable(cfg.hourly_file, ...
    'VariableNamingRule','preserve','TextType','string');
S = load(cfg.mat_file,'Q_net_raw','AreaByCell_m2');

required_hourly = { ...
    'UTC','ReferenceInterpolatedArea_km2', ...
    'CumulativePotentialIceReference_1e6m3', ...
    'CumulativePotentialIceConservative_1e6m3'};
assert(all(ismember(required_hourly,H.Properties.VariableNames)), ...
    'Hourly table does not contain all required columns.');
assert(isfield(S,'Q_net_raw') && isfield(S,'AreaByCell_m2'), ...
    'MAT file must contain Q_net_raw and AreaByCell_m2.');

if isdatetime(H.UTC)
    t = H.UTC;
    if isempty(t.TimeZone); t.TimeZone = 'UTC'; end
else
    t = datetime(string(H.UTC), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss', ...
        'TimeZone','UTC');
end

A_reference = double(H.ReferenceInterpolatedArea_km2);
V_reference = double(H.CumulativePotentialIceReference_1e6m3);
V_conservative = double(H.CumulativePotentialIceConservative_1e6m3);

assert(size(S.Q_net_raw,1) == numel(t), ...
    'MAT-file heat-flux time dimension does not match the hourly table.');
assert(size(S.Q_net_raw,2) == size(S.AreaByCell_m2,2), ...
    'ERA5-cell dimensions differ between Q_net_raw and AreaByCell_m2.');

%% 3. Continuous potential open-water forcing

% The fixed reference footprint is defined by the maximum expert-confirmed
% polynya area observed in each ERA5 cell during the study period. It is used
% only to summarize spatially varying ERA5 forcing; it is not the hourly
% polynya area used in the ice-production calculation.
fixed_cell_weights = max(S.AreaByCell_m2,[],1,'omitnan');
fixed_cell_weights(~isfinite(fixed_cell_weights)) = 0;
assert(sum(fixed_cell_weights) > 0, ...
    'No positive SAR-footprint weights were found in the MAT file.');

valid_q = isfinite(S.Q_net_raw);
q_numerator = sum( ...
    S.Q_net_raw .* fixed_cell_weights,2,'omitnan');
q_denominator = sum(double(valid_q).*fixed_cell_weights,2);
Q_openwater_hourly = q_numerator./q_denominator;
Q_openwater_hourly(q_denominator == 0) = NaN;

% A centered 24 h mean contains 11 preceding samples, the current sample,
% and 12 following samples.
Q_openwater_24h = movmean( ...
    Q_openwater_hourly,[11 12], ...
    'omitnan','Endpoints','shrink');

active = isfinite(A_reference) & A_reference > 0;
active_start = find(diff([false; active]) == 1);
active_end = find(diff([active; false]) == -1);
assert(numel(active_start) == numel(active_end), ...
    'Active interpolation windows could not be paired.');

salt_Mt_per_1e6m3 = ...
    cfg.rho_ice_kg_m3 * ...
    cfg.water_salinity_g_kg * ...
    cfg.salt_rejection_fraction / 1e6;

%% 4. Common time axis

t_start = dateshift(min(t),'start','month');
t_end = dateshift(max(t),'start','month')+calmonths(1);
month_ticks = (t_start:calmonths(1):t_end)';
english_months = [ ...
    "Jan","Feb","Mar","Apr","May","Jun", ...
    "Jul","Aug","Sep","Oct","Nov","Dec"];
month_labels = strings(size(month_ticks));
for k = 1:numel(month_ticks)
    month_name = english_months(month(month_ticks(k)));
    if k == 1 || month(month_ticks(k)) == 1
        month_labels(k) = sprintf('%s %d', ...
            month_name,year(month_ticks(k)));
    else
        month_labels(k) = month_name;
    end
end

%% 5. Figure 4

fig_main = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 cfg.main_width_cm+5 cfg.main_height_cm+5], ...
    'Renderer','painters');

tl = tiledlayout(fig_main,2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% Panel (a): potential open-water net heat flux

ax1 = nexttile(tl,1);
hold(ax1,'on');

q_all = [Q_openwater_hourly; Q_openwater_24h];
q_all = q_all(isfinite(q_all));
assert(~isempty(q_all),'No valid open-water heat-flux values were found.');
q_min = min([q_all; 0]);
q_max = max([q_all; 0]);
q_pad = max(25,0.06*max(q_max-q_min,1));
q_limits = [floor((q_min-q_pad)/50)*50 ...
            ceil((q_max+q_pad)/50)*50];
ylim(ax1,q_limits);

% Shade only SAR-supported interpolation windows. This identifies when the
% continuous potential forcing contributes to the area-weighted budget.
h_opening = gobjects(1);
for k = 1:numel(active_start)
    h_tmp = fill(ax1, ...
        [t(active_start(k)) t(active_end(k)) ...
         t(active_end(k)) t(active_start(k))], ...
        [q_limits(1) q_limits(1) q_limits(2) q_limits(2)], ...
        clr.opening, ...
        'FaceAlpha',cfg.opening_shade_alpha, ...
        'EdgeColor','none');
    if k == 1
        h_opening = h_tmp;
        h_opening.DisplayName = 'SAR-defined polynya periods';
    else
        h_tmp.HandleVisibility = 'off';
    end
end

h_q_hourly = plot(ax1,t,Q_openwater_hourly, ...
    '-', ...
    'Color',[0.66 0.67 0.69], ...
    'LineWidth',0.48, ...
    'DisplayName','Hourly net heat flux');
h_q_24h = plot(ax1,t,Q_openwater_24h, ...
    '-', ...
    'Color',clr.red, ...
    'LineWidth',cfg.line_width+0.15, ...
    'DisplayName','24 h mean net heat flux');

yline(ax1,0,'-', ...
    'Color',[0.38 0.38 0.38], ...
    'LineWidth',0.75, ...
    'HandleVisibility','off');
ylabel(ax1,'Potential open-water net heat loss (W m^{-2})', ...
    'FontSize',cfg.label_font_size);

if numel(active_start) > 0
    legend_handles = [h_q_24h h_q_hourly h_opening];
else
    legend_handles = [h_q_24h h_q_hourly];
end
lgd1 = legend(ax1,legend_handles, ...
    'Location','northwest', ...
    'NumColumns',3, ...
    'FontSize',cfg.legend_font_size, ...
    'Box','on');
lgd1.ItemTokenSize = [12 8];

%% Panel (b): cumulative potential ice production and salt rejection

ax2 = nexttile(tl,2);
hold(ax2,'on');
yyaxis(ax2,'left');

ice_span = max(V_reference-V_conservative,0);
h_envelope = area(ax2,t,[V_conservative ice_span], ...
    'LineStyle','none');
h_envelope(1).FaceColor = 'none';
h_envelope(1).FaceAlpha = 0;
h_envelope(1).HandleVisibility = 'off';
h_envelope(2).FaceColor = clr.orange;
h_envelope(2).FaceAlpha = cfg.band_alpha;
h_envelope(2).DisplayName = 'Temporal-area sensitivity range';

h_conservative = plot(ax2,t,V_conservative, ...
    '--', ...
    'Color',clr.grey, ...
    'LineWidth',cfg.sensitivity_line_width, ...
    'DisplayName','Conservative scenario');
h_reference = plot(ax2,t,V_reference, ...
    '-', ...
    'Color',clr.orange, ...
    'LineWidth',cfg.line_width+0.25, ...
    'DisplayName','Reference scenario');

ylabel(ax2,'Potential ice volume (10^6 m^3)', ...
    'FontSize',cfg.label_font_size);
ice_top = ceil((max(V_reference,[],'omitnan')+2)/5)*5;
ylim(ax2,[0 ice_top]);
ice_limits = ylim(ax2);
ax2.YAxis(1).Color = clr.orange;

yyaxis(ax2,'right');
ylim(ax2,ice_limits*salt_Mt_per_1e6m3);
ylabel(ax2,'Potential salt rejection (Mt)', ...
    'FontSize',cfg.label_font_size);
ax2.YAxis(2).Color = clr.dark;

yyaxis(ax2,'left');
lgd2 = legend(ax2,[h_reference h_conservative h_envelope(2)], ...
    'Location','northwest', ...
    'NumColumns',3, ...
    'FontSize',cfg.legend_font_size, ...
    'Box','on');
lgd2.ItemTokenSize = [12 8];

%% 6. Shared styling and export

linkaxes([ax1 ax2],'x');
axes_all = [ax1 ax2];
for k = 1:numel(axes_all)
    ax = axes_all(k);
    xlim(ax,[t_start t_end]);
    ax.XTick = month_ticks;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = clr.grid;
    ax.GridAlpha = 0.65;
    ax.Box = 'on';
    ax.Layer = 'top';
    ax.LineWidth = 0.75;
    ax.TickDir = 'out';
    ax.TickLength = [0.006 0.006];
    ax.FontName = cfg.font_name;
    ax.FontSize = cfg.axis_font_size;
    ax.XColor = clr.dark;
end

ax1.XTickLabel = [];
ax2.XTickLabel = cellstr(month_labels);
ax2.YAxis(1).Color = clr.orange;
ax2.YAxis(2).Color = clr.dark;

text(ax1,0.012,1.075,'(a)', ...
    'Units','normalized', ...
    'FontName',cfg.font_name, ...
    'FontWeight','bold', ...
    'FontSize',cfg.label_font_size, ...
    'VerticalAlignment','top', ...
    'Clipping','off');
text(ax2,0.012,1.075,'(b)', ...
    'Units','normalized', ...
    'FontName',cfg.font_name, ...
    'FontWeight','bold', ...
    'FontSize',cfg.label_font_size, ...
    'VerticalAlignment','top', ...
    'Clipping','off');

exportgraphics(fig_main,cfg.main_png, ...
    'Resolution',cfg.export_dpi, ...
    'BackgroundColor','white');
exportgraphics(fig_main,cfg.main_pdf, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

fprintf('\nFigure 4 heat-flux/ice-production figure exported.\n');
fprintf('PNG: %s\n',cfg.main_png);
fprintf('PDF: %s\n',cfg.main_pdf);
fprintf('Reference salt rejection fraction: %.2f\n', ...
    cfg.salt_rejection_fraction);
