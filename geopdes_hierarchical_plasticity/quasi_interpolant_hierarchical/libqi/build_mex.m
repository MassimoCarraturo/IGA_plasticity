function build_mex()
% BUILD_MEX  Compile both libqi MATLAB interfaces (qi_mex, qi_local_ls_mex).
%
%   Run this once from any folder; it cd's into the libqi root, invokes
%   `mex` with compiler-appropriate OpenMP flags and writes
%   matlab/qi_mex.<mexext>          (THB quasi-interpolant entry point)
%   matlab/qi_local_ls_mex.<mexext> (local LS B-spline projection entry)
%
%   On Windows with the Microsoft Visual C++ compiler, the OpenMP flag is
%   /openmp; on GCC/Clang it is -fopenmp. The script picks the right one.
%
%   Requirements: a C compiler set up via `mex -setup C`.

    here = fileparts(mfilename('fullpath'));
    old  = pwd;
    cleanupObj = onCleanup(@() cd(old));
    cd(here);

    cc = mex.getCompilerConfigurations('C', 'Selected');
    if isempty(cc)
        error('build_mex:nocompiler', ...
              'No C compiler is set up. Run "mex -setup C" first.');
    end

    isMSVC = contains(cc.Manufacturer, 'Microsoft', 'IgnoreCase', true);

    targets = struct();
    targets(1).name    = 'qi_mex';
    targets(1).entry   = 'matlab/qi_mex.c';
    targets(1).sources = {'src/qi.c', 'src/qi_linalg.c', 'src/qi_bspline.c'};

    targets(2).name    = 'qi_local_ls_mex';
    targets(2).entry   = 'matlab/qi_local_ls_mex.c';
    targets(2).sources = {'src/qi_localls.c', 'src/qi_linalg.c'};

    for k = 1:numel(targets)
        t = targets(k);
        if isMSVC
            compflags = 'COMPFLAGS=$COMPFLAGS /openmp /O2';
            mex('-R2018a', '-Iinclude', '-Isrc', compflags, ...
                t.entry, t.sources{:}, '-outdir', 'matlab');
        else
            cflags  = 'CFLAGS=$CFLAGS -fopenmp -O3 -std=c99';
            ldflags = 'LDFLAGS=$LDFLAGS -fopenmp';
            mex('-R2018a', '-Iinclude', '-Isrc', cflags, ldflags, ...
                t.entry, t.sources{:}, '-outdir', 'matlab');
        end
        fprintf('libqi MEX built -> %s\n', ...
                fullfile(here, 'matlab', [t.name, '.', mexext()]));
    end
end
