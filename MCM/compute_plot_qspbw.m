%% Plots of BWss
nnodes = 100;
% Read the position of the graphs
tol_ls=[0,10,20];
positions = dir('../KSP/graphs_qspbwss/star_roads/vx2/node_pos');  % List all graph position files
positions = positions(~ismember({positions.name}, {'.', '..'}));
positions = positions(contains({positions.name}, 'pos'));

accum_x = zeros(length(positions), nnodes);
accum_y = zeros(length(positions), nnodes);
accum_bw2 = zeros(length(positions), length(tol_ls), nnodes);

folders = dir('mcm/star_road_p_structure_vx2/');  % List all items in the 'mcm' directory
folders = folders([folders.isdir]);  % Keep only directories

% Remove '.' and '..' (current and parent directories)
folders = folders(~ismember({folders.name}, {'.', '..'}));

% Filter directories containing 'ring'
ring_folders = folders(contains({folders.name}, 'ring_road_'));
len_fold = length(ring_folders);


for file_i = 1: len_fold
    str_mat = sprintf(['./mcm/star_road_p_structure_vx2/%s/' ...
        'path_ls_mat_%d.mat'], ring_folders(file_i).name, 20);

    s = string(ring_folders(file_i).name);
    nums = double(extract(s, digitsPattern));  % -> [123 2]
    graph_num = nums(1);

    str_pos = sprintf(['../KSP/graphs_qspbwss/star_roads/vx2/' ...
        'node_pos/pos_ring%d.csv'], graph_num);



    T = readtable(str_pos, 'VariableNamingRule','preserve'); 

    % Access by column name:
    nodes = T.vertex1;
    x = T.x;
    y = T.y;
    
    %read cell array
    mat_struct = load(str_mat);
    fieldNames = fieldnames(mat_struct);

    if isempty(fieldNames)
        sprintf("entered %s\n", ring_folders(file_i).name);
        continue
    end
    PathLisArr = mat_struct.(fieldNames{1});

    if isempty(PathLisArr)
        continue
    end

    qsp_bw = compute_qsp_bw_from_file(str_mat, tol_ls , true, false, true, false);

    accum_x(graph_num+1,:) = x;
    accum_y(graph_num+1,:) = y;    
    accum_bw2(graph_num+1,:,:) = reshape(qsp_bw, [1, length(tol_ls), nnodes]);
    file_i
end

% write bwss to file

% Save the accumulated data to a .mat file for later analysis
save('./computed_bwss/accumulated_data_star.mat', 'accum_x', 'accum_y', 'accum_bw2');

%%

idx_all = ~all(accum_bw2==0, 3);

% accum_bw_noz = (accum_bw2(idx_all,:)-accum_bw(idx_all,:))/((nnodes-1)*(nnodes-2));
accum_bw_noz = accum_bw2(:, 2, idx_all);

x_nnz = accum_x(idx_all,:);
y_nnz = accum_y(idx_all,:);

accum_bw_noz = reshape(accum_bw_noz, [size(accum_bw_noz,1)*size(accum_bw_noz,3),1]);
x_nnz = reshape(x_nnz, [size(x_nnz,1)*size(x_nnz,2),1]);
y_nnz = reshape(y_nnz, [size(y_nnz,1)*size(y_nnz,2),1]);

% x, y, z are vectors (same length)

% --- 1) Get desired binning from binscatter (e.g., auto or fixed NumBins)
h = binscatter(x_nnz, y_nnz, 'NumBins', [40 40]);   % tweak bins to taste
xe = h.XBinEdges;                            % x-edges
ye = h.YBinEdges;                            % y-edges
delete(h);                                   % we’ll draw our own colored image

% --- 2) Bin indices and mean(z) per bin
[ix,~] = discretize(x_nnz, xe);
[iy,~] = discretize(y_nnz, ye);
m = ~isnan(ix) & ~isnan(iy) & ~isnan(accum_bw_noz);

nx = numel(xe)-1; ny = numel(ye)-1;
sumz = accumarray([ix(m) iy(m)], accum_bw_noz(m), [nx ny], @sum, 0);
cnt  = accumarray([ix(m) iy(m)], 1,     [nx ny], @sum, 0);
meanz = sumz ./ cnt;                       % mean z in each bin
% meanz(cnt==0) = NaN;                       % empty bins
meanz(isnan(meanz)) = 0;
% --- 3) Plot like a binscatter heatmap
xc = (xe(1:end-1) + xe(2:end))/2;          % bin centers
yc = (ye(1:end-1) + ye(2:end))/2;
min(accum_bw_noz)
imagesc(xc, yc, meanz.');                  % transpose to align x rows/y cols
set(gca, 'YDir','normal')                  % so y increases upward
axis equal tight
%axis([-1 1 -1 1]);

colormap(winter);

c = colorbar;%cool parula colormap(cool)
max(meanz,[],"all")
caxis([min(meanz,[],"all"), max(meanz,[],"all")]);
c.Label.String = '20% QSP-BW';
c.Label.FontSize = 22;
c.FontSize = 22;
ax = gca;
ax.FontSize = 16; 

%% Select tolerance slice (tol index 2)
tol_i = 3;

% bw: [nGraphs x nnodes]
bw = squeeze(accum_bw2(:, tol_i, :));

% Optional: keep only graphs that have at least one nonzero node
keep_graph = any(bw ~= 0, 2);

bw = bw(keep_graph, :);
xg = accum_x(keep_graph, :);
yg = accum_y(keep_graph, :);

% Keep only nonzero entries (node-level mask)
mask = (bw ~= 0) & ~isnan(bw) & ~isnan(xg) & ~isnan(yg);

% Vectorize
accum_bw_noz = bw(mask);
x_nnz = xg(mask);
y_nnz = yg(mask);

% Ensure column vectors
accum_bw_noz = accum_bw_noz(:);
x_nnz = x_nnz(:);
y_nnz = y_nnz(:);

% --- 1) Get desired binning from binscatter
h = binscatter(x_nnz, y_nnz, 'NumBins', [40 40]);
xe = h.XBinEdges;
ye = h.YBinEdges;
delete(h);

%--- 2) Bin indices and mean(z) per bin
[ix,~] = discretize(x_nnz, xe);
[iy,~] = discretize(y_nnz, ye);

m = ~isnan(ix) & ~isnan(iy) & ~isnan(accum_bw_noz);

nx = numel(xe)-1;
ny = numel(ye)-1;

sumz = accumarray([ix(m) iy(m)], accum_bw_noz(m), [nx ny], @sum, NaN);
cnt  = accumarray([ix(m) iy(m)], 1,              [nx ny], @sum, 0);

meanz = sumz ./ cnt;
meanz(cnt==0) = NaN;              % keep empty bins as NaN (better than 0)

% --- 3) Plot heatmap
xc = (xe(1:end-1) + xe(2:end))/2;
yc = (ye(1:end-1) + ye(2:end))/2;

imagesc(xc, yc, meanz.');
set(gca, 'YDir','normal')
axis equal tight

colormap(flipud(copper))
c = colorbar;

% Robust caxis (ignore NaNs)
caxis([min(meanz,[],'all','omitnan'), max(meanz,[],'all','omitnan')]);

cbar_label = sprintf('%d %% QSP-BW', tol_ls(tol_i));
c.Label.String = cbar_label;
c.Label.FontSize = 22;
c.FontSize = 22;

ax = gca;
ax.FontSize = 16;
