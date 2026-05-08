function coeffs = get_QI_coeffs_c(hspace, hmsh, data)
%GET_QI_COEFFS_C  Drop-in replacement for get_QI_coeffs.m (C/OpenMP backend).
%
%   coeffs = get_QI_coeffs_c(hspace, hmsh, data)
%
%   Inputs are GeoPDEs hierarchical structures.  The function builds the
%   plain-C structures expected by qi_mex and returns the same
%   (ndof x ncomp) coefficient matrix as the MATLAB reference.
%
%   Requirements:
%     - qi_mex MEX file built (see ../Makefile, target "mex")
%     - Path containing qi_mex on the MATLAB path

    ndof_per_level = double(hspace.ndof_per_level(:));

    % Active indices per level (MATLAB stores them already 1-based)
    active = cell(hspace.nlevels, 1);
    for lv = 1:hspace.nlevels
        if isfield(hspace, 'active') && numel(hspace.active) >= lv
            active{lv} = double(hspace.active{lv}(:));
        else
            active{lv} = (1:ndof_per_level(lv))';
        end
    end

    % Mesh: nel_dir per level
    nel_dir = cell(hmsh.nlevels, 1);
    for lv = 1:hmsh.nlevels
        nel_dir{lv} = double(hmsh.nel_dir{lv}(:)');
    end

    % Breaks per level (optional, but always passed when available)
    breaks = cell(hmsh.nlevels, 1);
    if isfield(hmsh, 'mesh_of_level')
        for lv = 1:hmsh.nlevels
            breaks{lv} = hmsh.mesh_of_level(lv).breaks;
        end
    elseif isfield(hmsh, 'breaks')
        for lv = 1:hmsh.nlevels
            breaks{lv} = hmsh.breaks{lv};
        end
    end

    % Build the lightweight structs expected by qi_mex
    hsp = struct();
    hsp.degree         = double(hspace.degree(:)');
    hsp.nlevels        = double(hspace.nlevels);
    hsp.ndof           = double(hspace.ndof);
    hsp.ndof_per_level = ndof_per_level;
    hsp.active         = active;

    hms = struct();
    hms.nel_dir        = nel_dir;
    if any(~cellfun(@isempty, breaks))
        hms.breaks     = breaks;
    end

    dat = struct();
    dat.x              = double(data.x(:));
    dat.y              = double(data.y(:));
    f                  = double(data.f);
    if size(f, 1) ~= numel(dat.x)
        f = f.';
    end
    dat.f              = f;

    coeffs = qi_mex(hsp, hms, dat);
end
