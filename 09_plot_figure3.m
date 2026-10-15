%% 09_plot_figure3.m
% Purpose: Reproduce manuscript Figure 3.
% Figure 3 shows variation in wind-ice coupling during 22-26 March 2026.
% Panel (a) compares median SAR-derived drift along the ESSiP opening normal
% with a 1.5% wind-scaling reference. Panel (b) shows the drift-to-wind ratio.
% Positive normal components indicate motion away from the grounded ice ridge.
% MATLAB: R2022b

clear; clc; close all;

%% 1. Paths and parameters
project_root = fileparts(mfilename('fullpath'));
processed_dir = fullfile(project_root,'processed');
drift_dir = fullfile(processed_dir,'ice_drift_event_20260322_20260326');

pair_file = fullfile(drift_dir, ...
    'SeaIceDrift_Event_PairSummary.csv');

wind_file = fullfile(drift_dir, ...
    'ERA5_Wind_By_SAR_Pair.csv');

area_file = fullfile(processed_dir, ...
    'SAR_constrained_heat_budget_ERA5_cells_expertQC', ...
    'SAR_Polynya_Observation_Area.csv');

output_dir = fullfile(project_root,'figures_output','Figure3_Wind_Ice_Coupling');

wind_scaling = 0.015;
seconds_per_day = 86400;
meters_per_km = 1000;
ms_to_kmd = seconds_per_day/meters_per_km;  % 86.4

if ~isfolder(output_dir)
    mkdir(output_dir);
end

required_files = {pair_file,wind_file,area_file};
for k = 1:numel(required_files)
    if ~isfile(required_files{k})
        error('Required file not found:\n%s',required_files{k});
    end
end

%% 2. Read and align the three data sources
P = readtable(pair_file,'VariableNamingRule','preserve');
W = readtable(wind_file,'VariableNamingRule','preserve');
A = readtable(area_file,'VariableNamingRule','preserve');

P.StartUTC = to_utc_datetime(P.StartUTC);
P.EndUTC   = to_utc_datetime(P.EndUTC);
W.StartUTC = to_utc_datetime(W.StartUTC);
W.EndUTC   = to_utc_datetime(W.EndUTC);
A.UTC      = to_utc_datetime(A.UTC);

% Keep the four consecutive SAR pairs from 22 to 26 March 2026.
event_start = datetime(2026,3,22,0,0,0,'TimeZone','UTC');
event_end   = datetime(2026,3,27,0,0,0,'TimeZone','UTC');
P = P(P.StartUTC>=event_start & P.EndUTC<event_end,:);
P = sortrows(P,'StartUTC');

if height(P)~=4
    error('Expected four SAR pairs during 22-26 March, but found %d.',height(P));
end

n = height(P);
wind_n_ms = nan(n,1);
area_start_km2 = nan(n,1);
area_end_km2 = nan(n,1);

time_tolerance = minutes(10);

for k = 1:n
    iw = nearest_time_pair(W.StartUTC,W.EndUTC, ...
        P.StartUTC(k),P.EndUTC(k),time_tolerance);
    wind_n_ms(k) = W.OpeningNormalWind_m_s(iw);

    ia0 = nearest_time_index(A.UTC,P.StartUTC(k),time_tolerance);
    ia1 = nearest_time_index(A.UTC,P.EndUTC(k),time_tolerance);
    area_start_km2(k) = A.TotalArea_km2(ia0);
    area_end_km2(k)   = A.TotalArea_km2(ia1);
end

%% 3. Derived diagnostics
ice_median_kmd = P.MedianOpeningComponent_km_day;
ice_q25_kmd = P.OpeningQ25_km_day;
ice_q75_kmd = P.OpeningQ75_km_day;

wind_n_kmd = wind_n_ms*ms_to_kmd;
reference_kmd = wind_scaling*wind_n_kmd;

coupling_pct = 100*ice_median_kmd./wind_n_kmd;
coupling_q25_pct = 100*ice_q25_kmd./wind_n_kmd;
coupling_q75_pct = 100*ice_q75_kmd./wind_n_kmd;

delta_area_km2 = area_end_km2-area_start_km2;

period_label = strings(n,1);
for k = 1:n
    period_label(k) = sprintf('%s-%s Mar', ...
        datestr(P.StartUTC(k),'dd'),datestr(P.EndUTC(k),'dd'));
end

DiagnosticTable = table( ...
    P.StartUTC,P.EndUTC,wind_n_ms,reference_kmd, ...
    ice_median_kmd,ice_q25_kmd,ice_q75_kmd, ...
    coupling_pct,coupling_q25_pct,coupling_q75_pct, ...
    area_start_km2,area_end_km2,delta_area_km2, ...
    'VariableNames',{ ...
    'StartUTC','EndUTC','OpeningNormalWind_m_s', ...
    'WindScaledReference_km_day', ...
    'ObservedMedianDrift_km_day','ObservedQ25_km_day', ...
    'ObservedQ75_km_day','CouplingRatio_percent', ...
    'CouplingQ25_percent','CouplingQ75_percent', ...
    'StartArea_km2','EndArea_km2','AreaChange_km2'});

disp('Figure 3 diagnostics:');
disp(DiagnosticTable);

writetable(DiagnosticTable,fullfile(output_dir, ...
    'Figure3_Wind_Ice_Coupling_Diagnostics.csv'));

%% 4. Figure style
font_name = 'Arial';
tick_font_size = 10.0;
label_font_size = 12.0;
annotation_font_size = 10;
panel_font_size = 11;

observed_color = [0.835 0.243 0.188];   % vermillion
reference_color = [0.000 0.447 0.698];  % blue
zero_color = [0.35 0.35 0.35];
grid_color = [0.84 0.84 0.84];
expand_bg = [0.945 0.972 0.992];
contract_bg = [0.992 0.958 0.950];

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 17.8 12], ...
    'Renderer','painters');

ax1 = axes(fig,'Position',[0.105 0.420 0.865 0.45]);
ax2 = axes(fig,'Position',[0.105 0.115 0.865 0.25]);

x = (1:n)';

%% 5. Panel (a): observed and wind-scaled drift
hold(ax1,'on');

all_y1 = [ice_q25_kmd;ice_q75_kmd;reference_kmd;0];
y1_min = min(all_y1,[],'omitnan');
y1_max = max(all_y1,[],'omitnan');
y1_span = max(y1_max-y1_min,1);
y1_lim = [y1_min-0.14*y1_span, y1_max+0.40*y1_span+1];

draw_interval_backgrounds(ax1,x,delta_area_km2,y1_lim, ...
    expand_bg,contract_bg);

yline(ax1,0,'-','Color',zero_color,'LineWidth',0.85, ...
    'HandleVisibility','off');

h_ref = plot(ax1,x,reference_kmd,'--d', ...
    'Color',reference_color, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',reference_color, ...
    'MarkerSize',6.0, ...
    'LineWidth',1.45);

h_obs = errorbar(ax1,x,ice_median_kmd, ...
    ice_median_kmd-ice_q25_kmd, ...
    ice_q75_kmd-ice_median_kmd, ...
    '-o', ...
    'Color',observed_color, ...
    'MarkerFaceColor',observed_color, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize',6.5, ...
    'LineWidth',1.55, ...
    'CapSize',7);

for k = 1:n
    text(ax1,x(k),ice_median_kmd(k), ...
        sprintf('  %.2f',ice_median_kmd(k)), ...
        'FontName',font_name,'FontSize',annotation_font_size, ...
        'Color',observed_color, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment',vertical_alignment(ice_median_kmd(k)));

    xn = (x(k)-0.5)/n;
    text(ax1,xn,0.965, ...
        sprintf('U_{wind,n} = %.2f m s^{-1}',wind_n_ms(k)), ...
        'Units','normalized','Interpreter','tex', ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'FontName',font_name,'FontSize',annotation_font_size, ...
        'Color',[0.18 0.18 0.18]);

    text(ax1,xn,0.865, ...
        sprintf('\\DeltaA = %+.1f km^2',delta_area_km2(k)), ...
        'Units','normalized','Interpreter','tex', ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'FontName',font_name,'FontSize',annotation_font_size, ...
        'FontWeight','bold', ...
        'Color',area_text_color(delta_area_km2(k), ...
        reference_color,observed_color));
end

text(ax1,0,1.12,'(a)', ...
    'Units','normalized','FontName',font_name, ...
    'FontSize',panel_font_size,'FontWeight','bold', ...
    'VerticalAlignment','top');



xlim(ax1,[0.5 n+0.5]);
ylim(ax1,y1_lim);
set(ax1,'XTick',x,'XTickLabel',[], ...
    'FontName',font_name,'FontSize',tick_font_size, 'FontWeight','bold',...
    'TickDir','out','LineWidth',0.75,'Box','on', ...
    'Layer','top','YGrid','on','XGrid','off', ...
    'GridColor',grid_color,'GridAlpha',0.55);
ylabel(ax1,'U_{ice,n} (km d^{-1}) ', ...
    'FontName',font_name,'FontSize',label_font_size,'FontWeight','bold');
legend(ax1,[h_obs h_ref], ...
    {'SAR-derived drift (median and IQR)', ...
    '1.5% wind-scaling reference'}, ...
    'Location','northoutside','Orientation','horizontal', ...
    'Box','off','FontName',font_name,'FontSize',tick_font_size);

%% 6. Panel (b): coupling ratio
hold(ax2,'on');

all_y2 = [coupling_q25_pct;coupling_q75_pct;1.5;0];
y2_min = min(all_y2,[],'omitnan');
y2_max = max(all_y2,[],'omitnan');
y2_span = max(y2_max-y2_min,0.5);
y2_lim = [y2_min-0.16*y2_span, y2_max+0.20*y2_span];

draw_interval_backgrounds(ax2,x,delta_area_km2,y2_lim, ...
    expand_bg,contract_bg);

yline(ax2,0,'-','Color',zero_color,'LineWidth',0.85, ...
    'HandleVisibility','off');

h_ratio_ref = yline(ax2,100*wind_scaling,'--', ...
    'Color',reference_color,'LineWidth',1.15, ...
    'Label','1.5% reference', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom', ...
    'FontName',font_name,'FontSize',annotation_font_size);

for k = 1:n
    plot(ax2,[x(k) x(k)],[0 coupling_pct(k)],'-', ...
        'Color',observed_color,'LineWidth',2.0, ...
        'HandleVisibility','off');
end

h_ratio = errorbar(ax2,x,coupling_pct, ...
    coupling_pct-coupling_q25_pct, ...
    coupling_q75_pct-coupling_pct, ...
    'o','LineStyle','none', ...
    'Color',observed_color, ...
    'MarkerFaceColor',observed_color, ...
    'MarkerEdgeColor','w', ...
    'MarkerSize',7.0, ...
    'LineWidth',1.45, ...
    'CapSize',7);

for k = 1:n
    offset = 0.055*y2_span;
    if coupling_pct(k)>=0
        y_text = coupling_pct(k)+offset;
        va = 'bottom';
    else
        y_text = coupling_pct(k)-offset;
        va = 'top';
    end

    text(ax2,x(k),y_text,sprintf('%+.2f%%',coupling_pct(k)), ...
        'HorizontalAlignment','center','VerticalAlignment',va, ...
        'FontName',font_name,'FontSize',annotation_font_size, ...
        'FontWeight','bold','Color',observed_color);
end

text(ax2,0,1.17,'(b)', ...
    'Units','normalized','FontName',font_name, ...
    'FontSize',panel_font_size,'FontWeight','bold', ...
    'VerticalAlignment','top');


% xlabel(ax2,'SAR image interval', ...
%     'FontName',font_name,'FontSize',label_font_size);

xlim(ax2,[0.5 n+0.5]);
ylim(ax2,y2_lim);
set(ax2,'XTick',x,'XTickLabel',cellstr(period_label), ...
    'FontName',font_name,'FontSize',tick_font_size,'FontWeight','bold', ...
    'TickDir','out','LineWidth',0.75,'Box','on', ...
    'Layer','top','YGrid','on','XGrid','off', ...
    'GridColor',grid_color,'GridAlpha',0.55);
ylabel(ax2,'U_{ice,n}/U_{wind,n} (%)', ...
    'FontName',font_name,'FontSize',label_font_size,'FontWeight','bold');
% Panel (b) is self-explanatory; no additional legend is required.
set(h_ratio_ref,'HandleVisibility','off');
set(h_ratio,'HandleVisibility','off');

%% 7. Export
set(findall(fig,'-property','FontName'),'FontName',font_name);

png_file = fullfile(output_dir, ...
    'Figure_3_Wind_Ice_Coupling.png');
pdf_file = fullfile(output_dir, ...
    'Figure_3_Wind_Ice_Coupling.pdf');

exportgraphics(fig,png_file, ...
    'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,pdf_file, ...
    'ContentType','vector','BackgroundColor','white');

fprintf('\nFigure 3 exported:\n%s\n%s\n',png_file,pdf_file);

%% Local functions
function t = to_utc_datetime(x)
    if isdatetime(x)
        t = x;
        t.TimeZone = 'UTC';
        return;
    end

    s = string(x);
    formats = { ...
        'yyyy-MM-dd''T''HH:mm:ss', ...
        'yyyy-MM-dd HH:mm:ss', ...
        'dd-MMM-yyyy HH:mm:ss'};

    for k = 1:numel(formats)
        try
            t = datetime(s,'InputFormat',formats{k},'TimeZone','UTC');
            if all(~isnat(t))
                return;
            end
        catch
        end
    end

    t = datetime(s,'TimeZone','UTC');
end

function idx = nearest_time_index(all_time,target_time,tolerance)
    [dt,idx] = min(abs(all_time-target_time));
    if isempty(idx) || isnan(seconds(dt)) || dt>tolerance
        error('No observation found within %g minutes of %s.', ...
            minutes(tolerance),string(target_time));
    end
end

function idx = nearest_time_pair(start_all,end_all,start_target,end_target,tolerance)
    mismatch = abs(start_all-start_target)+abs(end_all-end_target);
    [dt,idx] = min(mismatch);
    if isempty(idx) || isnan(seconds(dt)) || dt>2*tolerance
        error('No wind interval matches %s to %s.', ...
            string(start_target),string(end_target));
    end
end

function draw_interval_backgrounds(ax,x,delta_area,ylims,expand_color,contract_color)
    for k = 1:numel(x)
        if delta_area(k)>=0
            c = expand_color;
        else
            c = contract_color;
        end

        patch(ax, ...
            [x(k)-0.5 x(k)+0.5 x(k)+0.5 x(k)-0.5], ...
            [ylims(1) ylims(1) ylims(2) ylims(2)], ...
            c,'EdgeColor','none','HandleVisibility','off');
    end
end

function c = area_text_color(delta_area,expand_color,contract_color)
    if delta_area>=0
        c = expand_color;
    else
        c = contract_color;
    end
end

function va = vertical_alignment(value)
    if value>=0
        va = 'bottom';
    else
        va = 'top';
    end
end
