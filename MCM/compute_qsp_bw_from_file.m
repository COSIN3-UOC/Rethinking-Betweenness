% COMPUTE QASI SHORTEST PATH BETWEENNESS CENTRALITY FROM PATH LIST CELL
% ARRAY
%
%   qsp_bw = compute_qsp_bw_from_file(bw_fname, normalized, directed)
%
%   Computes node QSP-BW centrality given a matrix of qsp
%   lists loaded from file. Each cell {s,t} in path_carray contains the list
%   of nodes along a qsp from s to t (0-based node IDs).
%
%   INPUTS:
%       bw_fname   - string, filename containing path_carray (cell array)
%       normalized - logical, whether to normalize betweenness values
%       directed   - (optional) logical, true if graph is directed
%                    default = false
%
%   OUTPUT:
%       qsp_bw     - 1 x nnodes vector of betweenness centrality values
%
%   NOTES:
%       - Assumes path_carray{s,t} contains node IDs in 0-based indexing.
%       - Normalization factor differs for directed vs. undirected graphs.
%         Undirected: divide by (n-1)(n-2)
%         Directed:   divide by (n-1)(n-2)/2
%

function qsp_bw = compute_qsp_bw_from_file(bw_fname, tol, normalized, ...
                                            directed, padded, h5)

    % Handle optional input
    n_tols = length(tol);

    % Ensure tol is a column vector for consistent processing and sort it
   
    tol(tol>1) = tol(tol>1)/100;
    [tol, ~] = sort(tol);
 
    
    if h5 == false

        path_carray = load(bw_fname);
        fieldNames = fieldnames(path_carray);
    
        PathLisArr = path_carray.(fieldNames{1});

    else

        PathLisArr = h5paths_to_cell(bw_fname);

    end

    nnodes = length(PathLisArr);
    qsp_bw = zeros(n_tols, nnodes);

    for s = 1:nnodes
        s*100/nnodes;
        if directed
            start_t = 1;
        else
            start_t = s+1;
        end

        for t = start_t:nnodes

            if s == t
                continue
            end

            path_ls_st = PathLisArr{s,t};
            if padded
                n_paths = length(path_ls_st);
            else
                n_paths = numel(path_ls_st);
            end

            if n_paths > 0
                if padded
                    path_ls = path_ls_st;
                    max_lens = path_ls(1, end)*(1+tol);
                    % A is MxK: columns 1..K-1 are nodes, column K is path length (meters)
                    nodes = path_ls(:,1:end-1);                 % node columns only
                    nc = size(nodes,2);
                    
                    nz = nodes ~= 0;                       % mark nonzeros
                    % last nonzero column index per row (0 if none)
                    lastNZ = max( nz .* (1:nc), [], 2 );
                    
                    % mask for entries strictly AFTER the last nonzero in each row
                    mask1 = (1:nc) > lastNZ; % implicit expansion to Mxnc
                    % set only zeros in that tail to -1
                    nodes(mask1 & nodes==0) = -1;
                    
                    path_ls(:,1:end-1) = nodes; % put back; last column untouched
                else
                    max_lens = path_ls_st{1, 1}(1, end)*(1+tol);
                    % get the element in path_ls_st with the max length to
                    % create a padded array

                    % returns a numeric vector of lengths, then the maximum
                    lens = cellfun(@length, path_ls_st);
                    max_l_p = max(lens);

                    % row builder: take v(1:end-1), pad to length max_l_p with 0, append v(end) or NaN if empty
                    rowFromPath = @(v) [ ...
                        (numel(v)>=1) * v(1:min(numel(v)-1,max_l_p)), ...    % path
                        ones(1, max(0, max_l_p - max(0,numel(v)-1)))*-1, ...    % padding if needed (prealloc)
                        (numel(v)>1 && numel(v)-1 < max_l_p) * 0 + ...       % placeholder (keeps lengths consistent)
                        (numel(v)>=1) * v(end)...        % last element
                    ];
                    
                   rows = cellfun(rowFromPath, path_ls_st, 'UniformOutput', false);

                   % stack into matrix
                   path_ls = vertcat(rows{:});
                end
                % for node = 0:nnodes-1
                %     for max_len_i = 1:n_tols
                %         max_len = max_lens(max_len_i);
                %         if (node ~= (s-1) && node ~= (t-1))
                %             mask = path_ls(:, end) <= max_len; % rows to keep
                %             n_p_in_tol = sum(mask);
                %             counts_node = sum(path_ls(mask, 1:end-1) == node, 'all')...
                %             / n_p_in_tol;
                %             qsp_bw(max_len_i, node+1) = qsp_bw(max_len_i, node+1) + counts_node;
                %         end
                %     end
                % end
                nodesMat = path_ls(:,1:end-1);
                lensVec  = path_ls(:,end);
                isSorted = all(diff(lensVec) >= 0);
                n_p_in_tol = 0;
                for max_len_i = 1:n_tols
                    max_len = max_lens(max_len_i);
                    if isSorted
                        M = length(lensVec);
                        while n_p_in_tol < M && lensVec(n_p_in_tol+1) <= max_len
                            n_p_in_tol = n_p_in_tol + 1;
                        end
                        % admissibles = 1:k
                        if n_p_in_tol == 0, continue; end
                        vals = nodesMat(1:n_p_in_tol,:);          % MxK

                    else

                        mask = lensVec <= max_len;
                        n_p_in_tol = sum(mask);
                        if n_p_in_tol == 0, continue; end
                        vals = nodesMat(mask,:);          % MxK
                    end

                    vals = vals(vals >= 0);           % treu padding (-1) (i el 0 si és node real no el treguis)
                    cnt = accumarray(vals(:)+1, 1, [nnodes,1]);
                
                    cnt(s) = 0;       % perquè s és 1-based però node id = s-1 => index s
                    cnt(t) = 0;
                
                    qsp_bw(max_len_i,:) = qsp_bw(max_len_i,:) + (cnt.'/n_p_in_tol);
                end

            end
        end
    end

    % Normalization
    if normalized
        if directed
            qsp_bw = qsp_bw / ((nnodes-1) * (nnodes-2));
        else
            qsp_bw = 2 * qsp_bw / ((nnodes-1) * (nnodes-2));
        end
    end
end
