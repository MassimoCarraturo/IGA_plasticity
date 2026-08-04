function QI_coeff = getcoeff_localLS_Bspl_c(hspace, hmsh, data, f, lambda, weight)
%GETCOEFF_LOCALLS_BSPL_C  C/OpenMP drop-in for getcoeff_localLS_Bspl.
%
%   QI_coeff = getcoeff_localLS_Bspl_c(hspace, hmsh, data, f, lambda)
%   QI_coeff = getcoeff_localLS_Bspl_c(hspace, hmsh, data, f, lambda, weight)
%
%   Same algorithm as getcoeff_localLS_Bspl.m: per active hierarchical
%   B-spline, solve a small penalised local least-squares system formed by
%   collocation at scattered data points, with optional Tikhonov term on
%   the gradgrad mass matrix.  Identical signature; the per-DOF loop runs
%   in C with OpenMP parallelism.
%
%   Requirements:
%     - qi_local_ls_mex must be on the MATLAB path (build via libqi/build_mex.m)
%     - GeoPDEs must be available so we can call basisfun_multi and
%       op_gradgradu_gradgradv_tp once per level (these matrices depend only
%       on the level mesh, not on the per-DOF loop, so we precompute them
%       in MATLAB and hand them to libqi).

    if nargin < 6
        weight = [];
    end

    par_dim = numel(hspace.space_of_level(1).degree);
    deg     = double(hspace.space_of_level(1).degree(:)');

    if size(data, 2) ~= par_dim
        error('getcoeff_localLS_Bspl_c:dim', ...
              'data must have %d columns (par_dim)', par_dim);
    end
    [npts, ncomp] = size(f);
    if size(data, 1) ~= npts
        error('getcoeff_localLS_Bspl_c:size', ...
              'data and f must have the same number of rows');
    end

    nlevels = hspace.nlevels;

    % Per-level: pre-compute mass (gradgrad) and col (collocation) matrices.
    mass_cell  = cell(nlevels, 1);
    col_cell   = cell(nlevels, 1);
    ndof_dir   = zeros(nlevels, par_dim);
    nel_dir    = zeros(nlevels, par_dim);
    breaks     = cell(nlevels, 1);
    conn_first = cell(nlevels, 1);     % per-cell first-1D-basis index

    for lv = 1:nlevels
        sp_lv = hspace.space_of_level(lv);
        ms_lv = hmsh.mesh_of_level(lv);
        ndof_dir(lv, :) = double(sp_lv.ndof_dir(:)');
        nel_dir(lv, :)  = double(ms_lv.nel_dir(:)');
        breaks{lv}      = cell(par_dim, 1);
        conn_first{lv}  = cell(par_dim, 1);
        for d = 1:par_dim
            breaks{lv}{d} = double(ms_lv.breaks{d}(:)');
            % sp_univ(d).connectivity is (deg+1) x nel; row 1 gives the
            % 1-based index of the first 1D basis on each cell.
            conn_first{lv}{d} = double(sp_lv.sp_univ(d).connectivity(1, :));
        end
        if hspace.ndof_per_level(lv) > 0
            % qi_localls.c:325 reads mass only under `lambda != 0.0`; mirror that
            % predicate exactly, so a small nonzero qi_lambda keeps its Tikhonov term.
            if lambda ~= 0
                mass_cell{lv} = gradgrad_cached (sp_lv, ms_lv, lv, par_dim);
            else
                mass_cell{lv} = sparse(0, 0);
            end
            col_cell{lv}  = basisfun_multi(deg, data, sp_lv.knots);
        else
            % Empty placeholders (libqi never reads them when ndof_per_level=0)
            mass_cell{lv} = sparse(0, 0);
            col_cell{lv}  = sparse(0, npts);
        end
    end

    % Active indices: flatten level -> level cell array of 1-based linear indices
    %  -> single vector indexed by hierarchical DOF (k = 0..ndof-1).
    active_offset  = zeros(nlevels + 1, 1, 'int32');
    for lv = 1:nlevels
        active_offset(lv + 1) = active_offset(lv) + numel(hspace.active{lv});
    end
    if active_offset(end) ~= hspace.ndof
        error('getcoeff_localLS_Bspl_c:active_count', ...
              'sum of active counts (%d) does not match hspace.ndof (%d)', ...
              active_offset(end), hspace.ndof);
    end
    active_indices = zeros(hspace.ndof, 1);
    for lv = 1:nlevels
        rng = (active_offset(lv) + 1):active_offset(lv + 1);
        active_indices(rng) = double(hspace.active{lv}(:));
    end

    hsp = struct();
    hsp.par_dim         = par_dim;
    hsp.degree          = deg;
    hsp.nlevels         = nlevels;
    hsp.ndof            = double(hspace.ndof);
    hsp.ndof_per_level  = double(hspace.ndof_per_level(:));
    hsp.active_offset   = double(active_offset(:));
    hsp.active_indices  = active_indices;
    hsp.ndof_dir        = ndof_dir;
    hsp.nel_dir         = nel_dir;
    hsp.breaks          = breaks;
    hsp.conn_first      = conn_first;

    coords = double(data);
    fvals  = double(f);
    if isempty(weight)
        wvec = [];
    else
        wvec = double(weight(:));
    end

    QI_coeff = qi_local_ls_mex(hsp, coords, fvals, double(lambda), wvec, mass_cell, col_cell);
end

function M = gradgrad_cached (sp_lv, ms_lv, lv, par_dim)
% The gradgrad matrix depends only on the level space and mesh, both frozen once
% the level exists, but the projection is called up to 40 times per run. Cache it
% per level.
%
% The key must pin the geometry as well as the discretisation: the matrix is built
% from physical Hessians, so two runs with the same degree and subdivision on
% different geometries (cylinder vs sphere) would otherwise collide. Levels can
% also be rebuilt (hmsh_remove_empty_levels, num_bisections>0), so an index alone
% is not a key.
    persistent keys vals

    ctr  = num2cell (0.5 * ones (1, par_dim));
    gmap = reshape (double (ms_lv.map (ctr)), 1, []);
    key  = [double(sp_lv.degree(:)'), double(sp_lv.ndof_dir(:)'), ...
            double(ms_lv.nel_dir(:)'), cellfun(@numel, sp_lv.knots(:)'), ...
            cellfun(@(k) sum(k), sp_lv.knots(:)'), gmap];

    if numel(keys) >= lv && isequal (keys{lv}, key)
        M = vals{lv};
        return
    end

    M = op_gradgradu_gradgradv_tp (sp_lv, sp_lv, ms_lv);

    % Skip caching matrices big enough to matter against the solver's own
    % footprint; recomputing them is cheaper than paging.
    if nnz(M) * 12 < 2.5e8
        keys{lv} = key;
        vals{lv} = M;
    elseif numel(vals) >= lv
        keys{lv} = [];
        vals{lv} = [];
    end
end
