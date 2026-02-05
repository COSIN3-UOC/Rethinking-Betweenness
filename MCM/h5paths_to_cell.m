function C = h5paths_to_cell(fname, grp)
%H5PATHS_TO_CELL  Load path lists from HDF5 and build {i,j} cell-of-paths
%
%   C = H5PATHS_TO_CELL(fname)
%   C = H5PATHS_TO_CELL(fname, grp)
%
%   fname : path to the .h5 file
%   grp   : HDF5 group where 'paths', 'offsets', 'lengths', 'i', 'j', 'k' live
%           (default: '/paths')
%
%   Output:
%     C is an N×N cell array, where C{i,j} is a cell array of paths
%     (each path is a row vector of node IDs).

    if nargin < 2
        grp = '/paths';
    end
    info = h5info(fname, [grp '/paths']);
    info.Datatype

    % Read metadata and datasets
    N        = int32(h5readatt(fname, grp, 'N'));
    paths    = single(h5read(fname, [grp '/paths']));
    offsets  = int64(h5read(fname, [grp '/offsets']));
    lengths  = int32(h5read(fname, [grp '/lengths']));
    isrc     = int32(h5read(fname, [grp '/i']));
    jdst     = int32(h5read(fname, [grp '/j']));
    kidx     = int32(h5read(fname, [grp '/k']));  %#ok<NASGU>  % only used for sorting

    % Preallocate output: each C{i,j} will be a cell array of paths
    C = cell(N, N);

    % Build sorting order by (i, j, k) to preserve path order within each OD pair
    triplets = [isrc(:), jdst(:), kidx(:)];
    [~, order] = sortrows(triplets);   % order is a permutation of indices

    % Iterate in sorted order
    for n = 1:numel(order)
        idx = order(n);

        src = isrc(idx) + 1;   % convert 0-based to MATLAB 1-based
        dst = jdst(idx) + 1;

        start_pos = offsets(idx) + 1;           % HDF5 offsets are 0-based
        stop_pos  = offsets(idx) + int64(lengths(idx));

        try
            p = paths(start_pos:stop_pos);
        catch ME
            src
            dst
            start_pos
            stop_pos

            error('Stopping due to error: %s', ME.message);

        end

        if isempty(C{src, dst})
            C{src, dst} = {};
        end
        C{src, dst}{end+1, 1} = p.';
    end
end

