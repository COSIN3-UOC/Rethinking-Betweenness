%% PLOT RING ROAD PHASE TRANS VX1 (Figs 8 and 10)

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx1/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
D10_ALL = []; D20_ALL = []; 
Y0_ALL = []; Y10_ALL = []; Y20_ALL = [];

base_path = "./sim_results/all_graphs/penal/ring/vx1/"; % for fig 10
% base_path = "./sim_results/all_graphs/ring/vx1/"; %for fig 8

processing_rate = 1;

tol_ls = [0, 0.1, 0.2];
rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    % paths_path = paths_fname.folder+"/"+paths_fname.name+...
    % "/path_ls_mat_20.mat"; %for fig 8

    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/pen_path_ls_t_02_k_50_p_01.mat"; %for fig 10
    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
        normalized, directed, padded, false);
    
    % accept paths from i to j and j to i
    qsp_bw = 2*qsp_bw;
    nnodes = size(qsp_bw, 2);
    packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

    rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx1.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx1.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx1.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);

    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = sum(S0.DeltaNExp , 1);
    y10 = sum(S10.DeltaNExp, 1);
    y20 = sum(S20.DeltaNExp, 1);

    % --- differences vs tol0
    diff1 = y10 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
    diff2 = y20 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;

    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    % collect for optional aggregate lines
    D10_ALL = [D10_ALL; diff1];
    D20_ALL = [D20_ALL; diff2];

    Y0_ALL = [Y0_ALL; y0];
    Y10_ALL = [Y10_ALL; y10];
    Y20_ALL = [Y20_ALL; y20];
end

rho_c = rho_c./numel(graph_nums);
% Define colors

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;

% --- Optional: overlay median curves across graphs (thicker lines)
if ~isempty(D10_ALL)

    med10 = mean(D10_ALL, 1, 'omitnan');
    med20 = mean(D20_ALL, 1, 'omitnan');


    h10 = plot(possGenRates, med10, '-', 'Color', color10, 'LineWidth', 3, ...
   'HandleVisibility','on');
    
    h20 = plot(possGenRates, med20, '-', 'Color', color20, 'LineWidth', 3, ...
       'HandleVisibility','on');

    plot(possGenRates, zeros(length(possGenRates)), ...
        '--', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'HandleVisibility','off');

    errorbar(possGenRates, med10, std(D10_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color', color10, 'HandleVisibility','off');
    errorbar(possGenRates, med20, std(D20_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');
end



% ---- ONE legend only (main axes): eps=0,10,20,30 + black dot rho_c ----
ax = gca;

% dummy handle for legend
hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Dummy line for eps=0% (legend only)
h0 = plot(nan, nan, '-', 'Color', color0, 'LineWidth', 3);

% lgd = legend(ax, [h0 h10 h20 hRho], ...
%     {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', '$\varepsilon = 20\%$', '$\rho_c$'}, ...
%     'Location','southeast', 'Interpreter','latex');%for fig 8

lgd = legend(ax, [h0 h10 h20 hRho], ...
    {'$\hat{\varepsilon} = 0\%$', '$\hat{\varepsilon} = 10\%$', ...
    '$\hat{\varepsilon} = 20\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');%for fig 10

lgd.Position = [0.8 0.18, 0.10, 0.18];

ax = gca;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
xlabel('Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
ylabel('$\left \langle \eta - \eta_0 \right \rangle$','Interpreter','latex', 'FontSize', 20);

if exist('makePlotNice','file')==2
    scale = 1.2; %#ok<NASGU> 
    makePlotNice;
end

% --- inset (optional): show raw η

ax2 = axes('Position',[0.6 0.55 0.28 0.37]);
xlim(ax2,[0.03 3])
box on; hold(ax2,'on');
set(ax2, 'XScale', 'log');  % log-x for scatter

if ~isempty(graph_nums)
    g = graph_nums(1);

    med0 = mean(Y0_ALL, 1, 'omitnan');
    med10 = mean(Y10_ALL, 1, 'omitnan');
    med20 = mean(Y20_ALL, 1, 'omitnan');

    err0 = std(Y0_ALL./packetsGenPerMinute, 1, 'omitnan');
    err10 = std(Y10_ALL./packetsGenPerMinute, 1, 'omitnan');
    err20 = std(Y20_ALL./packetsGenPerMinute, 1, 'omitnan');


    plot(ax2, possGenRates, med0./ ...
        packetsGenPerMinute,...
        'Color', color0, 'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 0\%$');
    plot(ax2, possGenRates, med10./ ...
        packetsGenPerMinute,...
        'Color', color10,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 10\%$');
    plot(ax2, possGenRates, med20./ ...
        packetsGenPerMinute,...
        'Color', color20,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 20\%$');

    errorbar(ax2, possGenRates, med0./ ...
        packetsGenPerMinute , err0, 'LineStyle', ...
        'none','Color',color0, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med10./ ...
        packetsGenPerMinute, err10, 'LineStyle', ...
        'none','Color',color10, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med20./ ...
        packetsGenPerMinute, err20, 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');

    xline(ax2, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');

    xlabel(ax2,'Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
    ylabel(ax2,'$\langle \eta \rangle$','Interpreter','latex','FontSize',20);
end

% final sizing
width_inches  = 16; 
height_inches = 10;
set(gcf, 'Units','inches', 'Position',[0 0 width_inches height_inches]);
set(gcf, 'PaperUnits','inches', 'PaperSize',[width_inches height_inches]);

%% PLOT RING ROAD PHASE TRANS VX2

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx2/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
D10_ALL = []; D20_ALL = []; D30_ALL = [];
Y0_ALL = []; Y10_ALL = []; Y20_ALL = []; Y30_ALL = [];

firstLegend = true;  % only show three legend entries once
base_path = "./sim_results/all_graphs/ring/vx2/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2, 0.3];
% rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_30.mat";
    %"/pen_path_ls_t_03_k_50_p_01.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    % qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
    %     normalized, directed, padded, false);
    
    % accept paths from i to j and j to i
    % qsp_bw = 2*qsp_bw;
    nnodes = size(qsp_bw, 2);
    packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

    % rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    % rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx2.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx2.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx2.mat', g);
    file4 = sprintf(base_path+'tol30_%d_vx2.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);
    S30 = load(file4);

    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = sum(S0.DeltaNExp , 1);
    y10 = sum(S10.DeltaNExp, 1);
    y20 = sum(S20.DeltaNExp, 1);
    y30 = sum(S30.DeltaNExp, 1);

    % --- differences vs tol0
    diff1 = y10 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
    diff2 = y20 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
    diff3 = y30 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;

    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    % collect for optional aggregate lines
    D10_ALL = [D10_ALL; diff1];
    D20_ALL = [D20_ALL; diff2];
    D30_ALL = [D30_ALL; diff3];

    Y0_ALL = [Y0_ALL; y0];
    Y10_ALL = [Y10_ALL; y10];
    Y20_ALL = [Y20_ALL; y20];
    Y30_ALL = [Y30_ALL; y30];
end

% rho_c = rho_c./numel(graph_nums);
% Define colors

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;

% --- Optional: overlay median curves across graphs (thicker lines)
if ~isempty(D10_ALL)

    med10 = mean(D10_ALL, 1, 'omitnan');
    med20 = mean(D20_ALL, 1, 'omitnan');
    med30 = mean(D30_ALL, 1, 'omitnan');


    h10 = plot(possGenRates, med10, '-', 'Color', color10, 'LineWidth', 3, ...
   'HandleVisibility','on');
    
    h20 = plot(possGenRates, med20, '-', 'Color', color20, 'LineWidth', 3, ...
       'HandleVisibility','on');
    
    h30 = plot(possGenRates, med30, '-', 'Color', color30, 'LineWidth', 3, ...
       'HandleVisibility','on');


    plot(possGenRates, zeros(length(possGenRates)), ...
        '--', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'HandleVisibility','off');

    errorbar(possGenRates, med10, std(D10_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color', color10, 'HandleVisibility','off');
    errorbar(possGenRates, med20, std(D20_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');
    errorbar(possGenRates, med30, std(D30_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color',color30, 'HandleVisibility','off');
end



% ---- ONE legend only (main axes): eps=0,10,20,30 + black dot rho_c ----
ax = gca;

% dummy handle for legend
hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Dummy line for eps=0% (legend only)
h0 = plot(nan, nan, '-', 'Color', color0, 'LineWidth', 3);


% % Grab line handles from main axes by their DisplayName (you already set these)
% h0  = findobj(ax, 'Type','line', 'DisplayName','MC Simulation $\varepsilon = 0\%$');
% h10 = findobj(ax, 'Type','line', 'DisplayName','MC Simulation $\varepsilon = 10\%$');
% h20 = findobj(ax, 'Type','line', 'DisplayName','MC Simulation $\varepsilon = 20\%$');
% h30 = findobj(ax, 'Type','line', 'DisplayName','MC Simulation $\varepsilon = 30\%$');


lgd = legend(ax, [h0 h10 h20 h30 hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', '$\varepsilon = 20\%$', '$\varepsilon = 30\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.8 0.18, 0.10, 0.18];

ax = gca;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
xlabel('Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
ylabel('$\left \langle \eta - \eta_0 \right \rangle$','Interpreter','latex', 'FontSize', 20);

if exist('makePlotNice','file')==2
    scale = 1.2; %#ok<NASGU> 
    makePlotNice;
end

% --- inset (optional): show raw η

ax2 = axes('Position',[0.6 0.55 0.28 0.37]);
xlim(ax2,[0.03 3])
box on; hold(ax2,'on');
set(ax2, 'XScale', 'log');  % log-x for scatter

if ~isempty(graph_nums)
    g = graph_nums(1);

    med0 = mean(Y0_ALL, 1, 'omitnan');
    med10 = mean(Y10_ALL, 1, 'omitnan');
    med20 = mean(Y20_ALL, 1, 'omitnan');
    med30 = mean(Y30_ALL, 1, 'omitnan');

    err0 = std(Y0_ALL./packetsGenPerMinute, 1, 'omitnan');
    err10 = std(Y10_ALL./packetsGenPerMinute, 1, 'omitnan');
    err20 = std(Y20_ALL./packetsGenPerMinute, 1, 'omitnan');
    err30 = std(Y30_ALL./packetsGenPerMinute, 1, 'omitnan');

    % translucent markers with scatter
    % scatter(ax2, possGenRates, med0, 20, 'filled', ...
    %         'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor','none', ...
    %         'MarkerFaceAlpha', 0.5, 'DisplayName', '$\varepsilon=0\%$');
    % scatter(ax2, possGenRates, med10, 20, 'filled', ...
    %         'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor','none', ...
    %         'MarkerFaceAlpha', 0.5, 'DisplayName', '$\varepsilon=10\%$');
    % scatter(ax2, possGenRates, med20, 20, 'filled', ...
    %         'MarkerFaceColor', [1 1 0], 'MarkerEdgeColor','none', ...
    %         'MarkerFaceAlpha', 0.5, 'DisplayName', '$\varepsilon=20\%$');
    % scatter(ax2, possGenRates, med30, 36, 'filled', ...
    %         'MarkerFaceColor', [0.55 0.2 1], 'MarkerEdgeColor','none', ...
    %         'MarkerFaceAlpha', 0.5, 'DisplayName', '$\varepsilon=30\%$');

    plot(ax2, possGenRates, med0./ ...
        packetsGenPerMinute,...
        'Color', color0, 'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 0\%$');
    plot(ax2, possGenRates, med10./ ...
        packetsGenPerMinute,...
        'Color', color10,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 10\%$');
    plot(ax2, possGenRates, med20./ ...
        packetsGenPerMinute,...
        'Color', color20,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 20\%$');
    plot(ax2, possGenRates, med30./ ...
        packetsGenPerMinute,...
        'Color', color30,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 30\%$');

    errorbar(ax2, possGenRates, med0./ ...
        packetsGenPerMinute , err0, 'LineStyle', ...
        'none','Color',color0, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med10./ ...
        packetsGenPerMinute, err10, 'LineStyle', ...
        'none','Color',color10, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med20./ ...
        packetsGenPerMinute, err20, 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med30./ ...
        packetsGenPerMinute, err30, 'LineStyle', ...
        'none','Color',color30, 'HandleVisibility','off');

    xline(ax2, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(4), '--', 'Color', color30, 'LineWidth', 2, 'HandleVisibility','off');

    xlabel(ax2,'Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
    ylabel(ax2,'$\langle \eta \rangle$','Interpreter','latex','FontSize',20);
end

% final sizing
width_inches  = 16; 
height_inches = 10;
set(gcf, 'Units','inches', 'Position',[0 0 width_inches height_inches]);
set(gcf, 'PaperUnits','inches', 'PaperSize',[width_inches height_inches]);

%% PLOT CITY GRAPH (1 REALISATION) FIG 10

addpath('./mcm');

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;
% --- x-axis (generation rates) and packets/min
possGenRates = cat(2, linspace(0.0001,0.008,10), logspace(-2, log10(3), 100));      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on


base_path = "./sim_results/all_graphs/penal/bcn_full/";

file1 = sprintf(base_path+'tol0_vx2.mat');
file2 = sprintf(base_path+'tol10_vx2.mat');
file3 = sprintf(base_path+'tol20_vx2.mat');
file4 = sprintf(base_path+'tol30_vx2.mat');

file1_hires = sprintf(base_path+'tol0_vx2_hires.mat');
file2_hires = sprintf(base_path+'tol10_vx2_hires.mat');
file3_hires = sprintf(base_path+'tol20_vx2_hires.mat');
file4_hires = sprintf(base_path+'tol30_vx2_hires.mat');

need = {file1,file2,file3};%,file4};
if any(~cellfun(@(f) exist(f,'file')==2, need))
    warning('Missing file(s) for graph.');
end

S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
S10 = load(file2);
S20 = load(file3);
S30 = load(file4);

S0_hires  = load(file1_hires);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
S10_hires = load(file2_hires);
S20_hires = load(file3_hires);
S30_hires = load(file4_hires);

% --- sums across nodes → 1×Nrho
y0_end  = sum(S0.DeltaNExp , 1);
y10_end = sum(S10.DeltaNExp, 1);
y20_end = sum(S20.DeltaNExp, 1);
y30_end = sum(S30.DeltaNExp, 1);

y0_hires  = sum(S0_hires.DeltaNExp , 1);
y10_hires = sum(S10_hires.DeltaNExp, 1);
y20_hires = sum(S20_hires.DeltaNExp, 1);
y30_hires = sum(S30_hires.DeltaNExp, 1);

y0 = cat(2, y0_hires, y0_end);
y10 = cat(2, y10_hires, y10_end);
y20 = cat(2, y20_hires, y20_end);
y30 = cat(2, y30_hires, y30_end);


nnodes = length(S30.DeltaNExp);

packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

% --- differences vs tol0
diff1 = y10 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
diff2 = y20 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
diff3 = y30 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;

% --- plot on the same axes (translucent error bands via patch)
x = possGenRates;

% compute rho_c

paths_path = './mcm/other_graphs/bcn_amb_full/pen_path_ls_t_03_k_20_p_01.h5';

% compute qsp bwss to get the critical rhos
normalized=false;
directed = false;  % Set directed flag for graph processing
padded = false;    % Set padded flag for data processing
tol_ls = [0, 0.1, 0.2, 0.3];

% qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
%     normalized, directed, padded, true);
% qsp_bw = 2*qsp_bw;
% writematrix(qsp_bw, './sim_results/bcn_qspbwss_t_03_k_20_p_01.csv')
% processing_rate = 1;
% 
% rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
% rho_c = min(rhos_crit,[],2);
ax = gca;

% --- Optional: overlay median curves across graphs (thicker lines)
if ~isempty(diff1)

    h10 = plot(possGenRates, diff1, '-', 'Color', color10, 'LineWidth', ...
       3, 'HandleVisibility','on', 'DisplayName', ...
       'MC Simulation $\varepsilon = 10\%$');

    h20 = plot(possGenRates, diff2, '-', 'Color', color20, ...
        'LineWidth', 3, 'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 20\%$');

    h30 = plot(possGenRates, diff3, '-', 'Color', color30, 'LineWidth', 3, 'HandleVisibility', ...
        'on', 'DisplayName','MC Simulation $\varepsilon = 30\%$');

    plot(possGenRates, zeros(length(possGenRates)), ...
        '--', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'HandleVisibility','off');

    xline(ax, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
    xline(ax, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax, rho_c(4), '--', 'Color', color30, 'LineWidth', 2, 'HandleVisibility','off');

end
% dummy handle for legend
hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Dummy line for eps=0% (legend only)
h0 = plot(nan, nan, '-', 'Color', color0, 'LineWidth', 3);

lgd = legend(ax, [h0 h10 h20 h30 hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', '$\varepsilon = 20\%$', ...
    '$\varepsilon = 30\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.44 0.21, 0.10, 0.18];
% legend('Position',[0.68, 0.63],'Interpreter','latex');
%set(ax, 'XScale', 'log');  % log-x for scatter
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
xlabel('Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
ylabel('$\eta - \eta_0$','Interpreter','latex', 'FontSize', 20);
set(ax, 'XScale', 'log');  % log-x for scatter

if exist('makePlotNice','file')==2
    scale = 1.2; %#ok<NASGU>
    makePlotNice;
end

% --- inset (optional): show raw η

ax2 = axes('Position',[0.5807,0.1921,0.28 0.38]); %#ok<LAXES>
xlim(ax2,[0 3])
box on; hold(ax2,'on');
set(ax2, 'XScale', 'log');  % log-x for scatter


plot(ax2, possGenRates, y0./ ...
    packetsGenPerMinute,...
    'Color', color0, 'LineWidth', 3, ...
    'HandleVisibility','on', 'DisplayName', ...
    'MC Simulation $\varepsilon = 0\%$');
plot(ax2, possGenRates, y10./ ...
    packetsGenPerMinute,...
    'Color', color10,  'LineWidth', 3, ...
    'HandleVisibility','on', 'DisplayName', ...
    'MC Simulation $\varepsilon = 10\%$');
plot(ax2, possGenRates, y20./ ...
    packetsGenPerMinute,...
    'Color', color20,  'LineWidth', 3, ...
    'HandleVisibility','on', 'DisplayName', ...
    'MC Simulation $\varepsilon = 20\%$');
plot(ax2, possGenRates, y30./ ...
    packetsGenPerMinute,...
    'Color', color30,  'LineWidth', 3, ...
    'HandleVisibility','on', 'DisplayName', ...
    'MC Simulation $\varepsilon = 30\%$');

xline(ax2, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
xline(ax2, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax2, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax2, rho_c(4), '--', 'Color', color30, 'LineWidth', 2, 'HandleVisibility','off');


xlabel(ax2,'Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
ylabel(ax2,'$\eta$','Interpreter','latex', 'FontSize', 20);

% final sizing
width_inches  = 16; 
height_inches = 10;
set(gcf, 'Units','inches', 'Position',[0 0 width_inches height_inches]);
set(gcf, 'PaperUnits','inches', 'PaperSize',[width_inches height_inches]);


%% THROUGHPUT vs RHO vx1

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx1/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
thr0_ALL = []; thr10_ALL = []; thr20_ALL = [];
err0_ALL = []; err10_ALL = []; err20_ALL = [];

base_path = "./sim_results/all_graphs/ring/vx1/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2];
% rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_20.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    % qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
    %     normalized, directed, padded, false);
    % 
    % % accept paths from i to j and j to i
    % qsp_bw = 2*qsp_bw;
    % nnodes = size(qsp_bw, 2);
    % packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho
    % 
    % rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    % rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx1.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx1.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx1.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);


    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = S0.noutstat;
    y10 = S10.noutstat;
    y20 = S20.noutstat;
    
    err0 = S0.noutstat_std;
    err10 = S10.noutstat_std;
    err20 = S20.noutstat_std;

    y0 = y0(:,1);
    y10 = y10(:,1);
    y20 = y20(:,1);

    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    thr0_ALL = [thr0_ALL; y0.'];
    thr10_ALL = [thr10_ALL; y10.'];
    thr20_ALL = [thr20_ALL; y20.'];
   
    err0_ALL = [err0_ALL; err0.'];
    err10_ALL = [err10_ALL; err10.'];
    err20_ALL = [err20_ALL; err20.'];


end
% rho_c = rho_c./numel(graph_nums);

ax = gca;

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;

mean0 = mean(thr0_ALL, 1);
mean10 = mean(thr10_ALL, 1);
mean20 = mean(thr20_ALL, 1);

mean_err0 = sqrt(sum(err0_ALL.^2, 1))./numel(graph_nums);
mean_err10 = sqrt(sum(err10_ALL.^2, 1))./numel(graph_nums);
mean_err20 = sqrt(sum(err20_ALL.^2, 1))./numel(graph_nums);

% Plot the throughput data for different epsilon values
h0 = plot(ax, possGenRates, mean0, ...
    'Color', color0, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 0\%$');
h10 = plot(ax, possGenRates, mean10, ...
    'Color', color10, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 10\%$');
h20 = plot(ax, possGenRates, mean20, ...
    'Color', color20, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');
% 
% errorbar(ax, possGenRates, mean0 , mean_err0, 'LineStyle', ...
%         'none','Color',color0, 'HandleVisibility','off');
% errorbar(ax, possGenRates, mean10, mean_err10, 'LineStyle', ...
%     'none','Color',color10, 'HandleVisibility','off');
% errorbar(ax, possGenRates, mean20, mean_err20, 'LineStyle', ...
%     'none','Color',color20, 'HandleVisibility','off');

% --- 0 ---
fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean0 - mean_err0, fliplr(mean0 + mean_err0)], ...
    color0, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');
plot(ax, possGenRates, mean0, 'Color', color0, 'LineWidth', 1.5)

% --- 10 ---
fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean10 - mean_err10, fliplr(mean10 + mean_err10)], ...
    color10, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');
plot(ax, possGenRates, mean10, 'Color', color10, 'LineWidth', 1.5)

% --- 20 ---
fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean20 - mean_err20, fliplr(mean20 + mean_err20)], ...
    color20, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');
plot(ax, possGenRates, mean20, 'Color', color20, 'LineWidth', 1.5)


xline(ax, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');

xlabel('$\rho$','interpreter','latex', 'FontSize', 20);
ylabel('Throughput','interpreter','latex', 'FontSize', 20);

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
scale = 1.2;

hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Set the size of the figure (in inches)
width_inches = 16; % Width of the figure in inches
height_inches = 10; % Height of the figure in inches

% Convert inches to pixels (assuming 100 pixels per inch)
width_pixels = width_inches; % Width of the figure in pixels
height_pixels = height_inches; % Height of the figure in pixels

% Set the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width_inches, height_inches]); % [left, bottom, width, height]

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [width_inches, height_inches]);

set(ax, 'YScale', 'log');  % log-x for scatter
set(ax, 'XScale', 'log');  % log-x for scatter

lgd = legend(ax, [h0 h10 h20 hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', ...
    '$\varepsilon = 20\%$' '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.7 0.7, 0.10, 0.18];

makePlotNice;
%% THROUGHPUT vs RHO vx2

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx2/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
thr0_ALL = []; thr10_ALL = []; thr20_ALL = []; thr30_ALL = [];
err0_ALL = []; err10_ALL = []; err20_ALL = []; err30_ALL = [];

base_path = "./sim_results/all_graphs/ring/vx2/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2, 0.3];
rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_30.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
        normalized, directed, padded, false);

    % accept paths from i to j and j to i
    qsp_bw = 2*qsp_bw;
    nnodes = size(qsp_bw, 2);
    packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

    rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx2.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx2.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx2.mat', g);
    file4 = sprintf(base_path+'tol30_%d_vx2.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);
    S30 = load(file4);


    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = S0.noutstat;
    y10 = S10.noutstat;
    y20 = S20.noutstat;
    y30 = S30.noutstat;

    err0 = S0.noutstat_std;
    err10 = S10.noutstat_std;
    err20 = S20.noutstat_std;
    err30 = S30.noutstat_std;

    y0 = y0(:,1);
    y10 = y10(:,1);
    y20 = y20(:,1);
    y30 = y30(:,1);

    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    thr0_ALL = [thr0_ALL; y0.'];
    thr10_ALL = [thr10_ALL; y10.'];
    thr20_ALL = [thr20_ALL; y20.'];
    thr30_ALL = [thr30_ALL; y30.'];

    err0_ALL = [err0_ALL; err0.'];
    err10_ALL = [err10_ALL; err10.'];
    err20_ALL = [err20_ALL; err20.'];
    err30_ALL = [err30_ALL; err30.'];


end
rho_c = rho_c./numel(graph_nums);

ax = gca;

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;

mean0 = mean(thr0_ALL, 1);
mean10 = mean(thr10_ALL, 1);
mean20 = mean(thr20_ALL, 1);
mean30 = mean(thr30_ALL, 1);

mean_err0 = sqrt(sum(err0_ALL.^2, 1))./numel(graph_nums);
mean_err10 = sqrt(sum(err10_ALL.^2, 1))./numel(graph_nums);
mean_err20 = sqrt(sum(err20_ALL.^2, 1))./numel(graph_nums);
mean_err30 = sqrt(sum(err30_ALL.^2, 1))./numel(graph_nums);

% Plot the throughput data for different epsilon values
h0 = plot(ax, possGenRates, mean0, ...
    'Color', color0, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 0\%$');
h10 = plot(ax, possGenRates, mean10, ...
    'Color', color10, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 10\%$');
h20 = plot(ax, possGenRates, mean20, ...
    'Color', color20, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');
h30 = plot(ax, possGenRates, mean30, ...
    'Color', color30, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');


% --- 0 ---
fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean0 - mean_err0, fliplr(mean0 + mean_err0)], ...
    color0, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

% --- 10 ---
fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean10 - mean_err10, fliplr(mean10 + mean_err10)], ...
    color10, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

% --- 20 ---
fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean20 - mean_err20, fliplr(mean20 + mean_err20)], ...
    color20, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

fill(ax, ...
    [possGenRates, fliplr(possGenRates)], ...
    [mean30 - mean_err30, fliplr(mean30 + mean_err30)], ...
    color30, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

xline(ax, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(4), '--', 'Color', color30, 'LineWidth', 2, 'HandleVisibility','off');

xlabel('$\rho$','interpreter','latex', 'FontSize', 20);
ylabel('Throughput','interpreter','latex', 'FontSize', 20);

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
scale = 1.2;

hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Set the size of the figure (in inches)
width_inches = 16; % Width of the figure in inches
height_inches = 10; % Height of the figure in inches

% Convert inches to pixels (assuming 100 pixels per inch)
width_pixels = width_inches; % Width of the figure in pixels
height_pixels = height_inches; % Height of the figure in pixels

% Set the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width_inches, height_inches]); % [left, bottom, width, height]

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [width_inches, height_inches]);

set(ax, 'YScale', 'log');  % log-x for scatter
set(ax, 'XScale', 'log');  % log-x for scatter

lgd = legend(ax, [h0 h10 h20, h30, hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', ...
    '$\varepsilon = 20\%$', '$\varepsilon = 30\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.7 0.7, 0.10, 0.18];

makePlotNice;

%% PLOT FOR MACROSCOPIC FUNDAMENTAL DIAGRAM

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx2/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
thr0_ALL = []; thr10_ALL = []; thr20_ALL = []; thr30_ALL = [];
x0_ALL = []; x10_ALL = []; x20_ALL = []; x30_ALL = [];

err0_ALL = []; err10_ALL = []; err20_ALL = []; err30_ALL = [];

base_path = "./sim_results/all_graphs/ring/vx2/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2, 0.3];
% rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_30.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    % qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
    %     normalized, directed, padded, false);
    % 
    % % accept paths from i to j and j to i
    % qsp_bw = 2*qsp_bw;
    % nnodes = size(qsp_bw, 2);
    % packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho
    % 
    % rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    % rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_limq_vx2.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_limq_vx2.mat', g);
    file3 = sprintf(base_path+'tol20_%d_limq_vx2.mat', g);
    file4 = sprintf(base_path+'tol30_%d_limq_vx2.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);
    S30 = load(file4);


    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = S0.noutstat;
    y10 = S10.noutstat;
    y20 = S20.noutstat;
    y30 = S30.noutstat;

    x0 = S0.nstat;
    x10 = S10.nstat;
    x20 = S20.nstat;
    x30 = S30.nstat;

    err0 = S0.nstat_std;
    err10 = S10.nstat_std;
    err20 = S20.nstat_std;
    err30 = S30.nstat_std;

    y0 = y0(:,1);
    y10 = y10(:,1);
    y20 = y20(:,1);
    y30 = y30(:,1);

    % --- plot on the same axes (translucent error bands via patch)
   
    thr0_ALL = [thr0_ALL; y0.'];
    thr10_ALL = [thr10_ALL; y10.'];
    thr20_ALL = [thr20_ALL; y20.'];
    thr30_ALL = [thr30_ALL; y30.'];
    
    x0_ALL = [x0_ALL; x0.'];
    x10_ALL = [x10_ALL; x10.'];
    x20_ALL = [x20_ALL; x20.'];
    x30_ALL = [x30_ALL; x30.'];

    err0_ALL = [err0_ALL; err0.'];
    err10_ALL = [err10_ALL; err10.'];
    err20_ALL = [err20_ALL; err20.'];
    err30_ALL = [err30_ALL; err30.'];


end
% rho_c = rho_c./numel(graph_nums);

ax = gca;

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;

mean0 = mean(thr0_ALL, 1);
mean10 = mean(thr10_ALL, 1);
mean20 = mean(thr20_ALL, 1);
mean30 = mean(thr30_ALL, 1);

xmean0 = mean(x0_ALL, 1);
xmean10 = mean(x10_ALL, 1);
xmean20 = mean(x20_ALL, 1);
xmean30 = mean(x30_ALL, 1);

mean_err0 = sqrt(sum(err0_ALL.^2, 1))./numel(graph_nums);
mean_err10 = sqrt(sum(err10_ALL.^2, 1))./numel(graph_nums);
mean_err20 = sqrt(sum(err20_ALL.^2, 1))./numel(graph_nums);
mean_err30 = sqrt(sum(err30_ALL.^2, 1))./numel(graph_nums);

% Plot the throughput data for different epsilon values
h0 = plot(ax, xmean0, mean0, ...
    'Color', color0, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 0\%$');
h10 = plot(ax, xmean10, mean10, ...
    'Color', color10, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 10\%$');
h20 = plot(ax, xmean20, mean20, ...
    'Color', color20, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');
h30 = plot(ax, xmean30, mean30, ...
    'Color', color30, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');


% --- 0 ---
fill(ax, ...
    [xmean0, fliplr(xmean0)], ...
    [mean0 - mean_err0, fliplr(mean0 + mean_err0)], ...
    color0, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

% --- 10 ---
fill(ax, ...
    [xmean10, fliplr(xmean10)], ...
    [mean10 - mean_err10, fliplr(mean10 + mean_err10)], ...
    color10, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

% --- 20 ---
fill(ax, ...
    [xmean20, fliplr(xmean20)], ...
    [mean20 - mean_err20, fliplr(mean20 + mean_err20)], ...
    color20, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

fill(ax, ...
    [xmean30, fliplr(xmean30)], ...
    [mean30 - mean_err30, fliplr(mean30 + mean_err30)], ...
    color30, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'HandleVisibility','off');

% xline(ax, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
% xline(ax, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
% xline(ax, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');
% xline(ax, rho_c(4), '--', 'Color', color30, 'LineWidth', 2, 'HandleVisibility','off');

xlabel('Packets in the system','interpreter','latex', 'FontSize', 20);
ylabel('Throughput','interpreter','latex', 'FontSize', 20);

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
scale = 1.2;

% hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Set the size of the figure (in inches)
width_inches = 16; % Width of the figure in inches
height_inches = 10; % Height of the figure in inches

% Convert inches to pixels (assuming 100 pixels per inch)
width_pixels = width_inches; % Width of the figure in pixels
height_pixels = height_inches; % Height of the figure in pixels

% Set the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width_inches, height_inches]); % [left, bottom, width, height]

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [width_inches, height_inches]);

% set(ax, 'YScale', 'log');  % log-x for scatter
% set(ax, 'XScale', 'log');  % log-x for scatter

lgd = legend(ax, [h0 h10 h20, h30], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', ...
    '$\varepsilon = 20\%$', '$\varepsilon = 30\%$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.7 0.7, 0.10, 0.18];

makePlotNice;

%% PLOT LATTICE ROAD PHASE TRANS VX2 (APPENDIX FIGS)

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/lattice_road_p_structure_vx2/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
D10_ALL = []; D20_ALL = [];
Y0_ALL = []; Y10_ALL = []; Y20_ALL = [];

firstLegend = true;  % only show three legend entries once
base_path = "./sim_results/all_graphs/lattice/vx2/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2];
rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_20.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
        normalized, directed, padded, false);
    
    % accept paths from i to j and j to i
    qsp_bw = 2*qsp_bw;
    nnodes = size(qsp_bw, 2);
    packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

    rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx2.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx2.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx2.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);

    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = sum(S0.DeltaNExp , 1);
    y10 = sum(S10.DeltaNExp, 1);
    y20 = sum(S20.DeltaNExp, 1);

    % --- differences vs tol0
    diff1 = y10 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
    diff2 = y20 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;

    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    % collect for optional aggregate lines
    D10_ALL = [D10_ALL; diff1];
    D20_ALL = [D20_ALL; diff2];

    Y0_ALL = [Y0_ALL; y0];
    Y10_ALL = [Y10_ALL; y10];
    Y20_ALL = [Y20_ALL; y20];
end

rho_c = rho_c./numel(graph_nums);
% Define colors

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;

% --- Optional: overlay median curves across graphs (thicker lines)
if ~isempty(D10_ALL)

    med10 = mean(D10_ALL, 1, 'omitnan');
    med20 = mean(D20_ALL, 1, 'omitnan');


    h10 = plot(possGenRates, med10, '-', 'Color', color10, 'LineWidth', 3, ...
   'HandleVisibility','on');
    
    h20 = plot(possGenRates, med20, '-', 'Color', color20, 'LineWidth', 3, ...
       'HandleVisibility','on');


    plot(possGenRates, zeros(length(possGenRates)), ...
        '--', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'HandleVisibility','off');

    errorbar(possGenRates, med10, std(D10_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color', color10, 'HandleVisibility','off');
    errorbar(possGenRates, med20, std(D20_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');

end

% ---- ONE legend only (main axes): eps=0,10,20,30 + black dot rho_c ----
ax = gca;

% dummy handle for legend
hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Dummy line for eps=0% (legend only)
h0 = plot(nan, nan, '-', 'Color', color0, 'LineWidth', 3);

lgd = legend(ax, [h0 h10 h20 hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', '$\varepsilon = 20\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.8 0.18, 0.10, 0.18];

ax = gca;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
xlabel('Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
ylabel('$\left \langle \eta - \eta_0 \right \rangle$','Interpreter','latex', 'FontSize', 20);

if exist('makePlotNice','file')==2
    scale = 1.2; %#ok<NASGU> 
    makePlotNice;
end

% --- inset (optional): show raw η

ax2 = axes('Position',[0.6 0.55 0.28 0.37]);
xlim(ax2,[0.03 3])
box on; hold(ax2,'on');
set(ax2, 'XScale', 'log');  % log-x for scatter

if ~isempty(graph_nums)
    g = graph_nums(1);

    med0 = mean(Y0_ALL, 1, 'omitnan');
    med10 = mean(Y10_ALL, 1, 'omitnan');
    med20 = mean(Y20_ALL, 1, 'omitnan');

    err0 = std(Y0_ALL./packetsGenPerMinute, 1, 'omitnan');
    err10 = std(Y10_ALL./packetsGenPerMinute, 1, 'omitnan');
    err20 = std(Y20_ALL./packetsGenPerMinute, 1, 'omitnan');

    plot(ax2, possGenRates, med0./ ...
        packetsGenPerMinute,...
        'Color', color0, 'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 0\%$');
    plot(ax2, possGenRates, med10./ ...
        packetsGenPerMinute,...
        'Color', color10,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 10\%$');
    plot(ax2, possGenRates, med20./ ...
        packetsGenPerMinute,...
        'Color', color20,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 20\%$');

    errorbar(ax2, possGenRates, med0./ ...
        packetsGenPerMinute , err0, 'LineStyle', ...
        'none','Color',color0, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med10./ ...
        packetsGenPerMinute, err10, 'LineStyle', ...
        'none','Color',color10, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med20./ ...
        packetsGenPerMinute, err20, 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');

    xline(ax2, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');

    xlabel(ax2,'Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
    ylabel(ax2,'$\langle \eta \rangle$','Interpreter','latex','FontSize',20);
end

% final sizing
width_inches  = 16; 
height_inches = 10;
set(gcf, 'Units','inches', 'Position',[0 0 width_inches height_inches]);
set(gcf, 'PaperUnits','inches', 'PaperSize',[width_inches height_inches]);

%% PLOT STAR ROAD PHASE TRANS VX2 (APPENDIX FIGS)

% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/star_road_p_structure_vx2/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
D10_ALL = []; D20_ALL = [];
Y0_ALL = []; Y10_ALL = []; Y20_ALL = [];

firstLegend = true;  % only show three legend entries once
base_path = "./sim_results/all_graphs/star/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2];
% rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    gi/numel(graph_nums)

    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_20.mat";

    if isempty(fieldnames(load(paths_path)))
        continue% struct has no fields
    end    % compute qsp bwss to get the critical rhos

    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    % qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
    %     normalized, directed, padded, false);
    % 
    % % accept paths from i to j and j to i
    % qsp_bw = 2*qsp_bw;
    % nnodes = size(qsp_bw, 2);
    % packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho
    % 
    % rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    % rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx2.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx2.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx2.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);

    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = sum(S0.DeltaNExp , 1);
    y10 = sum(S10.DeltaNExp, 1);
    y20 = sum(S20.DeltaNExp, 1);

    % --- differences vs tol0
    diff1 = y10 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;
    diff2 = y20 ./ packetsGenPerMinute - y0 ./ packetsGenPerMinute;

    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    % collect for optional aggregate lines
    D10_ALL = [D10_ALL; diff1];
    D20_ALL = [D20_ALL; diff2];

    Y0_ALL = [Y0_ALL; y0];
    Y10_ALL = [Y10_ALL; y10];
    Y20_ALL = [Y20_ALL; y20];
end

% rho_c = rho_c./numel(graph_nums);
% Define colors

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;

% --- Optional: overlay median curves across graphs (thicker lines)
if ~isempty(D10_ALL)

    med10 = mean(D10_ALL, 1, 'omitnan');
    med20 = mean(D20_ALL, 1, 'omitnan');


    h10 = plot(possGenRates, med10, '-', 'Color', color10, 'LineWidth', 3, ...
   'HandleVisibility','on');
    
    h20 = plot(possGenRates, med20, '-', 'Color', color20, 'LineWidth', 3, ...
       'HandleVisibility','on');


    plot(possGenRates, zeros(length(possGenRates)), ...
        '--', 'Color', [0 0 0 0.2], 'LineWidth', 2, 'HandleVisibility','off');

    errorbar(possGenRates, med10, std(D10_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color', color10, 'HandleVisibility','off');
    errorbar(possGenRates, med20, std(D20_ALL, 1, 'omitnan'), 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');

end

% ---- ONE legend only (main axes): eps=0,10,20,30 + black dot rho_c ----
ax = gca;

% dummy handle for legend
hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Dummy line for eps=0% (legend only)
h0 = plot(nan, nan, '-', 'Color', color0, 'LineWidth', 3);

lgd = legend(ax, [h0 h10 h20 hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', '$\varepsilon = 20\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.8 0.18, 0.10, 0.18];

ax = gca;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
xlabel('Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
ylabel('$\left \langle \eta - \eta_0 \right \rangle$','Interpreter','latex', 'FontSize', 20);

if exist('makePlotNice','file')==2
    scale = 1.2; %#ok<NASGU> 
    makePlotNice;
end

% --- inset (optional): show raw η

ax2 = axes('Position',[0.6 0.55 0.28 0.37]);
xlim(ax2,[0.03 3])
box on; hold(ax2,'on');
set(ax2, 'XScale', 'log');  % log-x for scatter

if ~isempty(graph_nums)
    g = graph_nums(1);

    med0 = mean(Y0_ALL, 1, 'omitnan');
    med10 = mean(Y10_ALL, 1, 'omitnan');
    med20 = mean(Y20_ALL, 1, 'omitnan');

    err0 = std(Y0_ALL./packetsGenPerMinute, 1, 'omitnan');
    err10 = std(Y10_ALL./packetsGenPerMinute, 1, 'omitnan');
    err20 = std(Y20_ALL./packetsGenPerMinute, 1, 'omitnan');

    plot(ax2, possGenRates, med0./ ...
        packetsGenPerMinute,...
        'Color', color0, 'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 0\%$');
    plot(ax2, possGenRates, med10./ ...
        packetsGenPerMinute,...
        'Color', color10,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 10\%$');
    plot(ax2, possGenRates, med20./ ...
        packetsGenPerMinute,...
        'Color', color20,  'LineWidth', 3, ...
        'HandleVisibility','on', 'DisplayName', ...
        'MC Simulation $\varepsilon = 20\%$');

    errorbar(ax2, possGenRates, med0./ ...
        packetsGenPerMinute , err0, 'LineStyle', ...
        'none','Color',color0, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med10./ ...
        packetsGenPerMinute, err10, 'LineStyle', ...
        'none','Color',color10, 'HandleVisibility','off');
    errorbar(ax2, possGenRates, med20./ ...
        packetsGenPerMinute, err20, 'LineStyle', ...
        'none','Color',color20, 'HandleVisibility','off');

    xline(ax2, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
    xline(ax2, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');

    xlabel(ax2,'Generation rate $\rho$','Interpreter','latex', 'FontSize', 20);
    ylabel(ax2,'$\langle \eta \rangle$','Interpreter','latex','FontSize',20);
end

% final sizing
width_inches  = 16; 
height_inches = 10;
set(gcf, 'Units','inches', 'Position',[0 0 width_inches height_inches]);
set(gcf, 'PaperUnits','inches', 'PaperSize',[width_inches height_inches]);


%% PLOT INDIVIDUAL PHASE TRANSITIONS
folders = dir('mcm/ring_road_p_structure_vx2/');  % List all items in the 'mcm' directory
folders = folders([folders.isdir]);  % Keep only directories

% Remove '.' and '..' (current and parent directories)
folders = folders(~ismember({folders.name}, {'.', '..'}));

% Filter directories containing 'ring'
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = length(ring_folders);

graph_nums = ones(len_fold)*-1;
for file_i=1:len_fold
    s = string(ring_folders(file_i).name);
    nums = double(extract(s, digitsPattern));  % -> [123 2]
    graph_num = nums(1);
    graph_nums(file_i) = graph_num;
end

addpath('./mcm');
numNodes = 100;
possGenRates = logspace(-2,log10(3),100);
for rhosIndex = 1:length(possGenRates) 
    % setting node generation rates
    genRatePerMinuteGates = zeros(numNodes,1);
    genRatePerMinuteGates(:) = possGenRates(rhosIndex);       
    packetsGenPerMinute(rhosIndex) = sum(genRatePerMinuteGates);
end

for file_i=1:length(graph_nums)
    graph_n = graph_nums(file_i);
    if graph_n == -1
        continue
    end

    file1 = sprintf("./sim_results/all_graphs/ring/vx2/tol0_%d_vx2.mat", ...
        graph_n);
    file2 = sprintf("./sim_results/all_graphs/ring/vx2/tol10_%d_vx2.mat", ...
    graph_n);
    file3 = sprintf("./sim_results/all_graphs/ring/vx2/tol20_%d_vx2.mat", ...
    graph_n);
    file4 = sprintf("./sim_results/all_graphs/ring/vx2/tol30_%d_vx2.mat", ...
    graph_n);

    DeltaNExp0 = load(file1);
    DeltaNExp10 = load(file2);
    DeltaNExp20 = load(file3);
    DeltaNExp30 = load(file4);
    
    diff1 = sum(DeltaNExp10.DeltaNExp) ./ (packetsGenPerMinute) -...
    sum(DeltaNExp0.DeltaNExp) ./ (packetsGenPerMinute);
    
    diff2 = sum(DeltaNExp20.DeltaNExp) ./ (packetsGenPerMinute) -...
    sum(DeltaNExp0.DeltaNExp) ./ (packetsGenPerMinute);
    
    diff3 = sum(DeltaNExp30.DeltaNExp) ./ (packetsGenPerMinute) -...
    sum(DeltaNExp0.DeltaNExp) ./ (packetsGenPerMinute);
    
    fig = figure ('Position', [50 50 650 500]);
    % semilogx(possGenRates,sum(DeltaNTeor) ./ (packetsGenPerMinute),'b-','LineWidth',2); hold on;
    % semilogx(possGenRates,sum(DeltaNExp0.DeltaNExp) ./ (packetsGenPerMinute),'bo','LineWidth',2); hold on;
    % semilogx(possGenRates,sum(DeltaNExp10.DeltaNExp) ./ (packetsGenPerMinute),'ro','LineWidth',2); hold on;
    % semilogx(possGenRates,sum(DeltaNExp20.DeltaNExp) ./ (packetsGenPerMinute),'yo','LineWidth',2); hold on;
    % semilogx(possGenRates,sum(DeltaNExp30.DeltaNExp) ./ (packetsGenPerMinute),'mo','LineWidth',2); hold on;
    % semilogx(possGenRates,sum(DeltaNExp40.DeltaNExp) ./ (packetsGenPerMinute),'co','LineWidth',2); hold on;
    % semilogx(possGenRates,sum(DeltaNExp50.DeltaNExp) ./ (packetsGenPerMinute),'go','LineWidth',2); hold on;
    
    
    semilogx(possGenRates,diff1,'ro','LineWidth',2); hold on;
    semilogx(possGenRates,diff2,'yo','LineWidth',2); hold on;
    semilogx(possGenRates,diff3,'mo','LineWidth',2); hold on;
    
    % semilogx(possGenRates,sum(DeltaNExp30.DeltaNExp) ./ (packetsGenPerMinute),'mo','LineWidth',2); hold on;
    
    err1 = sqrt(sum(DeltaNExp0.DeltaNExp_std.^2) + ...
        sum(DeltaNExp10.DeltaNExp_std.^2))./ (packetsGenPerMinute);
    err2 = sqrt(sum(DeltaNExp0.DeltaNExp_std.^2) + ...
        sum(DeltaNExp20.DeltaNExp_std.^2))./ (packetsGenPerMinute);
    err3 = sqrt(sum(DeltaNExp0.DeltaNExp_std.^2) + ...
        sum(DeltaNExp30.DeltaNExp_std.^2))./ (packetsGenPerMinute);
    
    errorbar(possGenRates,diff1, ...
        err1, 'vertical', 'LineStyle', 'none', 'color', 'r');
    errorbar(possGenRates,diff2, ...
        err2, 'vertical', 'LineStyle', 'none', 'color', 'y');
    errorbar(possGenRates,diff3, ...
        err3, 'vertical', 'LineStyle', 'none', 'color', 'm');

    xlabel('Generation rate $\rho$','interpreter','latex');
    ylabel('$\eta - \eta_0$','interpreter','latex');
    % legend('MC Simulation $\varepsilon = 0\%$',...
    %        'MC Simulation $\varepsilon = 10\%$',...
    %        'MC Simulation $\varepsilon = 20\%$',...
    %        'Location','northwest','interpreter','latex');
    %        %'Theoretical prediction',...
    %        % 'MC Simulation $\varepsilon = 40\%$',...
    %        % 'Critical injection rate $\rho_c$',...
    %        % '$\rho_c$ for $\varepsilon = 10\%$',...
    %        % '$\rho_c$ for $\varepsilon = 20\%$',...
    %        % '$\rho_c$ for $\varepsilon = 30\%$',...
    %        % '$\rho_c$ for $\varepsilon = 40\%$',...
    
    legend('MC Simulation $\varepsilon = 10\%$',...
           'MC Simulation $\varepsilon = 20\%$',...
           'MC Simulation $\varepsilon = 30\%$',...
           'Location','northwest','interpreter','latex');
    % axis([min(possGenRates) max(possGenRates) -0.07 0.12]);
    axis([min(possGenRates) max(possGenRates) -0.05 0.09]);
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
    yaxisproperties= get(gca, 'YAxis');
    yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
    scale = 1.2;
    makePlotNice;
    
    % crit_range = (possGenRates > rho_c0 - rho_c0/5) & (possGenRates < rho_c0+rho_c0/5);
    crit_range = (possGenRates > 0.06) & (possGenRates < 0.12);
    
    % create smaller axes in top right, and plot on it
    ax2 = axes('Position',[0.68 .6 .2 .3]);
    box on
    % y1 = sum(DeltaNTeor);
    y2 = sum(DeltaNExp0.DeltaNExp);
    y3 = sum(DeltaNExp10.DeltaNExp);
    y4 = sum(DeltaNExp20.DeltaNExp);
    % y5 = sum(DeltaNExp30.DeltaNExp);
    % y6 = sum(DeltaNExp40.DeltaNExp);
    % y7 = sum(DeltaNExp50.DeltaNExp);
    
    % semilogx(ax2, possGenRates(crit_range),y1(crit_range) ./ (packetsGenPerMinute(crit_range)),'b-','LineWidth',2); hold on;
    semilogx(ax2, possGenRates(crit_range),y2(crit_range) ./ ...
        (packetsGenPerMinute(crit_range)),'bo','LineWidth',2); hold on;
    semilogx(ax2, possGenRates(crit_range),y3(crit_range) ./ ...
        (packetsGenPerMinute(crit_range)),'ro','LineWidth',2); hold on;
    semilogx(ax2, possGenRates(crit_range),y4(crit_range) ./ ...
        (packetsGenPerMinute(crit_range)),'yo','LineWidth',2); hold on;

    xlabel('Generation rate $\rho$','interpreter','latex');
    ylabel('Order parameter $\eta$','interpreter','latex');
    % Set the size of the figure (in inches)
    width_inches = 16; % Width of the figure in inches
    height_inches = 10; % Height of the figure in inches
    legend('MC Simulation $\varepsilon = 0\%$',...
           'MC Simulation $\varepsilon = 10\%$',...
           'MC Simulation $\varepsilon = 20\%$',...
           'Location','northwest','interpreter','latex');
end
% Convert inches to pixels (assuming 100 pixels per inch)
width_pixels = width_inches; % Width of the figure in pixels
height_pixels = height_inches; % Height of the figure in pixels

% Set the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width_inches, height_inches]); % [left, bottom, width, height]

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [width_inches, height_inches]);

% Save the figure to PDF
% print('ring_tols.pdf', '-dpdf'); % '-r300' sets the resolution to 300 DPI (optional)


%% PLOT TOTAL DELIVERED PACKETS vx1
% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx1/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
thr0_ALL = []; thr10_ALL = []; thr20_ALL = [];

base_path = "./sim_results/all_graphs/ring/vx1/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2, 0.3];
rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_20.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
        normalized, directed, padded, false);

    % accept paths from i to j and j to i
    qsp_bw = 2*qsp_bw;
    nnodes = size(qsp_bw, 2);
    packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

    rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx1.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx1.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx1.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);

    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = S0.DelivExp;
    y10 = S10.DelivExp;
    y20 = S20.DelivExp;


    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
 
    thr0_ALL = [thr0_ALL; y0(:).'];
    thr10_ALL = [thr10_ALL; y10(:).'];
    thr20_ALL = [thr20_ALL; y20(:).'];
end
rho_c = rho_c./numel(graph_nums);

ax = gca;

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;

% Plot the throughput data for different epsilon values
h0 = plot(ax, possGenRates, mean(thr0_ALL, 1), ...
    'Color', color0, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 0\%$');
h10 = plot(ax, possGenRates, mean(thr10_ALL, 1), ...
    'Color', color10, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 10\%$');

h20 = plot(ax, possGenRates, mean(thr20_ALL, 1), ...
    'Color', color20, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');


xline(ax, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');


xlabel('Generation rate $\rho$','interpreter','latex', 'FontSize', 20);
ylabel('Total delivered packets','interpreter','latex', 'FontSize', 20);

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
scale = 1.2;
% crit_range = (possGenRates > rho_c0 - rho_c0/5) & (possGenRates < rho_c0+rho_c0/5);

% create smaller axes in top right, and plot on it
% ax2 = axes('Position',[.37 .2 .21 .21]);
% box on
% 
% scatter(ax2, rho_c0,0,150,'o', 'filled', 'b'); hold on;
% scatter(ax2, rho_c10,0,150, 'o', 'filled', 'r'); hold on;
% scatter(ax2, rho_c20,0,150, 'o', 'filled', 'y'); hold on;
% scatter(ax2, rho_c30,0,150, 'o', 'filled', 'm'); hold on;
% scatter(ax2, rho_c40,0,150, 'o', 'filled', 'c'); hold on;

hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Set the size of the figure (in inches)
width_inches = 16; % Width of the figure in inches
height_inches = 10; % Height of the figure in inches

% Convert inches to pixels (assuming 100 pixels per inch)
width_pixels = width_inches; % Width of the figure in pixels
height_pixels = height_inches; % Height of the figure in pixels

% Set the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width_inches, height_inches]); % [left, bottom, width, height]

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [width_inches, height_inches]);

set(ax, 'YScale', 'log');  % log-x for scatter
% set(ax, 'XScale', 'log');  % log-x for scatter

lgd = legend(ax, [h0 h10 h20  hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', '$\varepsilon = 20\%$', '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.7 0.7, 0.10, 0.18];

makePlotNice;

% Save the figure to PDF
% print('throughput_ring_slog.pdf', '-dpdf'); % '-r300' sets the resolution to 300 DPI (optional)
%% PLOT TOTAL DELIVERED PACKETS vx2
% --- find ring folders and extract numeric IDs robustly
folders = dir('mcm/ring_road_p_structure_vx2/');
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.','..'}));
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = numel(ring_folders);

graph_nums = -1 * ones(len_fold,1);
for k = 1:len_fold
    s = string(ring_folders(k).name);
    nums = double(extract(s, digitsPattern));   % may be empty
    if ~isempty(nums)
        graph_nums(k) = nums(1);
    end
end
graph_nums = graph_nums(graph_nums~=-1);   % drop missing

addpath('./mcm');

% --- x-axis (generation rates) and packets/min
possGenRates = logspace(-2, log10(3), 100);      % 1×Nrho

% Prepare the figure once
figure('Position',[50 50 900 600]); 
hold on; box on

% (Optional) collect diffs to overlay median/thick lines later
thr0_ALL = []; thr10_ALL = []; thr20_ALL = []; thr30_ALL = [];

base_path = "./sim_results/all_graphs/ring/vx2/";

processing_rate = 1;

tol_ls = [0, 0.1, 0.2, 0.3];
rho_c = zeros(length(tol_ls),1);

for gi = 1:numel(graph_nums)
    g = graph_nums(gi);
    
    pat = "_" + string(g);            % pattern: underscore followed by g
    
    % example when field is 'name'
    names = {ring_folders.name};           % 1×N cell of char (or string scalars)
    paths_fname = ring_folders(endsWith(string(names), pat));
    paths_path = paths_fname.folder+"/"+paths_fname.name+...
    "/path_ls_mat_30.mat";

    % compute qsp bwss to get the critical rhos
    normalized=false;
    directed = false;  % Set directed flag for graph processing
    padded = true;    % Set padded flag for data processing

    qsp_bw = compute_qsp_bw_from_file(paths_path, tol_ls, ...
        normalized, directed, padded, false);

    % accept paths from i to j and j to i
    qsp_bw = 2*qsp_bw;
    nnodes = size(qsp_bw, 2);
    packetsGenPerMinute = nnodes * possGenRates;   % 1×Nrho

    rhos_crit = processing_rate*(nnodes-1)./(qsp_bw+2*(nnodes-1));
    rho_c = rho_c+min(rhos_crit,[],2);

    file1 = sprintf(base_path+'tol0_%d_vx2.mat',  g);
    file2 = sprintf(base_path+'tol10_%d_vx2.mat', g);
    file3 = sprintf(base_path+'tol20_%d_vx2.mat', g);
    file4 = sprintf(base_path+'tol30_%d_vx2.mat', g);

    need = {file1,file2,file3};%,file4};
    if any(~cellfun(@(f) exist(f,'file')==2, need))
        warning('Missing file(s) for graph %d. Skipping.', g);
        continue
    end

    S0  = load(file1);   % expects: DeltaNExp (nodes×rho), DeltaNExp_std (nodes×rho)
    S10 = load(file2);
    S20 = load(file3);
    S30 = load(file4);


    % --- sums across nodes → 1×Nrho
    % If your arrays are (rho×nodes), use sum(...,2).' consistently instead.
    y0  = S0.DelivExp;
    y10 = S10.DelivExp;
    y20 = S20.DelivExp;
    y30 = S30.DelivExp;
    
    y0 = y0(:,1);
    y10 = y10(:,1);
    y20 = y20(:,1);
    y30 = y30(:,1);


    % --- plot on the same axes (translucent error bands via patch)
    x = possGenRates;
   
    thr0_ALL = [thr0_ALL; y0.'];
    thr10_ALL = [thr10_ALL; y10.'];
    thr20_ALL = [thr20_ALL; y20.'];
    thr30_ALL = [thr30_ALL; y30.'];

end
rho_c = rho_c./numel(graph_nums);

ax = gca;

color0 = [16 11 29]./255;
color10 = [186 89 151]./255;
color20 = [232 135 112]./255;
color30 = [250 201 62]./255;

% Plot the throughput data for different epsilon values
h0 = plot(ax, possGenRates, mean(thr0_ALL, 1), ...
    'Color', color0, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 0\%$');
h10 = plot(ax, possGenRates, mean(thr10_ALL, 1), ...
    'Color', color10, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 10\%$');
h20 = plot(ax, possGenRates, mean(thr20_ALL, 1), ...
    'Color', color20, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 20\%$');
h30 = plot(ax, possGenRates, mean(thr30_ALL, 1), ...
    'Color', color30, 'LineWidth', 3, 'HandleVisibility', 'on', ...
    'DisplayName', 'MC Simulation $\varepsilon = 30\%$');


xline(ax, rho_c(1), '--', 'Color', color0,  'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(2), '--', 'Color', color10, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(3), '--', 'Color', color20, 'LineWidth', 2, 'HandleVisibility','off');
xline(ax, rho_c(4), '--', 'Color', color30, 'LineWidth', 2, 'HandleVisibility','off');

xlabel('Generation rate $\rho$','interpreter','latex', 'FontSize', 20);
ylabel('Total delivered packets','interpreter','latex', 'FontSize', 20);

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
scale = 1.2;

hRho = plot(nan, nan, '--k', 'LineWidth', 2);

% Set the size of the figure (in inches)
width_inches = 16; % Width of the figure in inches
height_inches = 10; % Height of the figure in inches

% Convert inches to pixels (assuming 100 pixels per inch)
width_pixels = width_inches; % Width of the figure in pixels
height_pixels = height_inches; % Height of the figure in pixels

% Set the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width_inches, height_inches]); % [left, bottom, width, height]

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [width_inches, height_inches]);

set(ax, 'YScale', 'log');  % log-x for scatter
% set(ax, 'XScale', 'log');  % log-x for scatter

lgd = legend(ax, [h0 h10 h20 h30 hRho], ...
    {'$\varepsilon = 0\%$', '$\varepsilon = 10\%$', ...
    '$\varepsilon = 20\%$', '$\varepsilon = 30\%$' '$\rho_c$'}, ...
    'Location','southeast', 'Interpreter','latex');

lgd.Position = [0.7 0.7, 0.10, 0.18];

makePlotNice;

% Save the figure to PDF
% print('throughput_ring_slog.pdf', '-dpdf'); % '-r300' sets the resolution to 300 DPI (optional)