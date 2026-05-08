function build_mex()
% BUILD_MEX  Compile the qi_mex MATLAB interface to libqi.
%
%   Run this once from any folder; it cd's into the libqi root, invokes
%   `mex` with compiler-appropriate OpenMP flags and writes
%   matlab/qi_mex.mexw64 (or .mexa64 / .mexmaci64).
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

    if isMSVC
        compflags = 'COMPFLAGS=$COMPFLAGS /openmp /O2';
        mex('-R2018a', '-Iinclude', '-Isrc', compflags, ...
            'matlab/qi_mex.c', 'src/qi.c', 'src/qi_linalg.c', 'src/qi_bspline.c', ...
            '-outdir', 'matlab');
    else
        cflags  = 'CFLAGS=$CFLAGS -fopenmp -O3 -std=c99';
        ldflags = 'LDFLAGS=$LDFLAGS -fopenmp';
        mex('-R2018a', '-Iinclude', '-Isrc', cflags, ldflags, ...
            'matlab/qi_mex.c', 'src/qi.c', 'src/qi_linalg.c', 'src/qi_bspline.c', ...
            '-outdir', 'matlab');
    end

    fprintf('libqi MEX built successfully -> %s\n', ...
            fullfile(here, 'matlab', ['qi_mex.', mexext()]));
end
