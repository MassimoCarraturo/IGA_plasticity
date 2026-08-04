% STUDY_FRONT_METRICS_SPHERE  Front-focused metrics for the projection
% operators on the internally pressurised eighth-sphere (exact Hill reference).
%
% Spherical symmetry makes the von Mises condition sigma_t - sigma_r = sigma_y,
% so the yield constant here is PLAIN sigma_y (no 2/sqrt(3) factor, unlike the
% plane-strain cylinder).  The exact solution has
%   plastic (r<c): sigma_t = 2 s_y (1/2 - ln(c/r) - (1-c^3/b^3)/3), d/dr = +2 s_y/r
%   elastic (r>c): sigma_t = 2 s_y c^3/(3 b^3) (b^3/(2 r^3) + 1),   d/dr = -s_y c^3/r^4
% so the hoop-stress slope jumps by [d sigma_t/dr] = 3 s_y / c at the front.
%
% The mesh is driven by the GEOMETRIC front estimator ('plastic_front_sphere'),
% which is solution-independent, so all operators share an identical
% front-refined mesh.  Metric definitions: front_kink_metrics.m.
%
% QI-h (num_bisections = 1) is omitted: bisecting the projection mesh in 3-D
% multiplies its DOFs by eight and is not affordable here.
%
% Output: results/front_metrics_sphere/

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
results_dir = fullfile(here, 'results', 'front_metrics_sphere');
if ~exist(results_dir,'dir'); mkdir(results_dir); end

addpath(genpath(fullfile(project_root,'nurbs-1.4.3','nurbs-1.4.3','inst')));
addpath(genpath(fullfile(project_root,'geopdes-3.2.2','geopdes','inst')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','hierarchical_classes')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','adaptivity_iga')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','initialize')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','examples','plasticity')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','quasi_interpolant_hierarchical')));
addpath(fullfile(project_root,'geopdes_hierarchical_plasticity','quasi_interpolant_hierarchical','libqi','matlab'));
cd(here);
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');

%% physical data (as in study_sphere_qi_variants.m)
E=210000; nu=0.49; sigma_y=240; a=100; b=200; nload=5;
P = 332*0.99;                       % ~ 98.8 % of the limit pressure

pd.geo_name='geo_eighth_sphere.txt';
pd.nmnn_sides=[2]; pd.drchlt_sides=[]; pd.press_sides=[1];
pd.symm_sides=[3 4 5]; pd.slider_sides=[6];
pd.penalty_slider=@(x,y,z) 1e8*ones(size(x));
pd.yield_stress=@(x,y,z) sigma_y*ones(size(x));
pd.kappa_lame=@(x,y,z) E/(3*(1-2*nu))*ones(size(x));
pd.mu_lame=@(x,y,z) E/(2*(1+nu))*ones(size(x));
pd.f=@(x,y,z) zeros(3,size(x,1),size(x,2),size(x,3));
pd.g=@(x,y,z,ind) zeros(3,size(x,1),size(x,2),size(x,3));
pd.h=@(x,y,z,ind) zeros(3,size(x,1),size(x,2),size(x,3));
pd.p=@(x,y,z) P*ones(size(x));
% required by the geometric front estimator
pd.R_i=a; pd.R_o=b; pd.s_y=sigma_y; pd.Pmax=P;

c = fzero(@(x) -P + 2*sigma_y*log(x/a) + (2/3)*sigma_y*(1-x.^3/b^3), [a b]);
fprintf('front c = %.4f mm ; exact slope jump 3 s_y/c = %.4f MPa/mm\n', c, 3*sigma_y/c);

%% discretization: geometric front marking -> identical mesh for all operators
p=4;
md0.degree=[p p p]; md0.regularity=[p-1 p-1 p-1];
md0.nsub_coarse=[5 5 5]; md0.nsub_refine=[2 2 2];
md0.nquad=[p+1 p+1 p+1]; md0.space_type='standard'; md0.truncated=1; md0.nload=nload;
md0.newton_tol=1e-8; md0.newton_tol_abs=1e-10; md0.newton_iter_max=100;
ad0.flag='elements'; ad0.estimator='plastic_front_sphere'; ad0.C0_est=1.0;
ad0.mark_param=0.8; ad0.mark_param_coarsening=0.1; ad0.mark_strategy='MS';
ad0.max_ndof=1e7; ad0.max_nel=1e6; ad0.num_max_iter=3; ad0.tol=1e-10;   % must not bind: max_level is the control
ad0.adm_strategy='admissible'; ad0.coarsening_flag='any'; ad0.adm=p;

dband=12; djump=9; nsamp=3001;
method_tags={'L2','QI','QI_C0','Bezier','Bezier_C0','DLSQ','DLSQ_W'};   % QI_graded dropped from the study
max_levels=[2 3];

% Resume support.  A level-3 operator costs hours at p = 4, so a machine sleep
% or a stopped run must not discard finished work: metrics.mat is written after
% every operator (below) and reloaded here, and completed entries are skipped.
R=struct();
ckpt = fullfile(results_dir,'metrics.mat');
if exist(ckpt,'file')
    S_ck = load(ckpt,'R');
    if isfield(S_ck,'R'); R = S_ck.R; end
    fprintf('resuming: %d entries already in %s\n', numel(fieldnames(R)), ckpt);
end
for il=1:numel(max_levels)
for im=1:numel(method_tags)
    md=md0; ad=ad0; ad.max_level=max_levels(il);
    key = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
    % Skip only when BOTH the metric entry and its profile dump exist, so a rerun
    % after the profile output was added regenerates the missing profiles.
    pfile = fullfile(results_dir,'profiles',sprintf('%s_lev%d.dat',method_tags{im},max_levels(il)));
    if isfield(R, key) && exist(pfile,'file')
        fprintf('\n==== %s (max_level %d) : cached, skipped ====\n', method_tags{im}, max_levels(il));
        continue
    end
    switch method_tags{im}
        case 'L2';        md.type_projection='L2';
        case 'QI';        md.type_projection='QI';
        case 'QI_C0';     md.type_projection='QI_C0';
        case 'QI_graded'; md.type_projection='QI_graded'; md.graded_transition_level=2;
        case 'Bezier';    md.type_projection='BEZIER';
        case 'Bezier_C0'; md.type_projection='BEZIER_C0';
        case 'DLSQ';      md.type_projection='DLSQ';
        case 'DLSQ_W';    md.type_projection='DLSQ_W';
    end
    fprintf('\n==== %s (max_level %d) ====\n', method_tags{im}, max_levels(il));
    t0=tic;
    % One operator failing must not discard the other eight: record and continue.
    try
        [geo,chm,chsp,chs,~,~,csig,~,~]=adaptivity_J2_plasticity(pd,md,ad);
    catch ME
        fprintf('  SKIPPED %s (max_level %d): %s\n', method_tags{im}, max_levels(il), ME.message);
        continue
    end
    wall=toc(t0);

    hs=chs{end}; sig=csig{end};
    uu=linspace(0,1,nsamp); v0=0.5; w0=0.5;         % interior radial ray
    [~,F]=sp_eval(sig(:,1),hs,geo,{uu,v0,w0});
    X=reshape(F(1,:,:,:),[],1); Yc=reshape(F(2,:,:,:),[],1); Z=reshape(F(3,:,:,:),[],1);
    r=sqrt(X.^2+Yc.^2+Z.^2);
    er=[X./r, Yc./r, Z./r];                          % unit radial direction

    S=cell(1,6); G=cell(1,6);
    for k=1:6
        S{k}=reshape(sp_eval(sig(:,k),hs,geo,{uu,v0,w0}),[],1);
        g=sp_eval(sig(:,k),hs,geo,{uu,v0,w0},'gradient');
        G{k}=reshape(g(1,:,:,:),[],1).*er(:,1) + reshape(g(2,:,:,:),[],1).*er(:,2) ...
           + reshape(g(3,:,:,:),[],1).*er(:,3);      % d/dr of component k
    end
    % sigma_r = e_r' sigma e_r ; sigma_t = (tr sigma - sigma_r)/2  (sigma_t=sigma_phi)
    q = @(A) A{1}.*er(:,1).^2 + A{2}.*er(:,2).^2 + A{3}.*er(:,3).^2 ...
           + 2*A{4}.*er(:,1).*er(:,2) + 2*A{5}.*er(:,1).*er(:,3) + 2*A{6}.*er(:,2).*er(:,3);
    sr_h  = q(S);          tr_h  = S{1}+S{2}+S{3};
    dsr_h = q(G);          dtr_h = G{1}+G{2}+G{3};   % e_r is constant along the ray
    st_h  = 0.5*(tr_h  - sr_h);
    dst_h = 0.5*(dtr_h - dsr_h);
    vm_h  = sqrt(0.5*((S{1}-S{2}).^2+(S{2}-S{3}).^2+(S{3}-S{1}).^2)+3*(S{4}.^2+S{5}.^2+S{6}.^2));

    [st_e,dst_e] = hill_sph_t(r,b,c,sigma_y);
    M = front_kink_metrics(r, st_h, dst_h, st_e, dst_e, vm_h, sigma_y, c, dband, djump);
    % Profile dump for the stress-kink figures: sigma_t and its analytic
    % counterpart along the sampling ray, subsampled to keep the files small.
    pdir = fullfile(results_dir,'profiles');
    if ~exist(pdir,'dir'); mkdir(pdir); end
    ss = max(1, round(numel(r)/800));
    idx = 1:ss:numel(r);
    Pm = [r(idx), st_h(idx), st_e(idx), sr_h(idx)];
    fid_p = fopen(fullfile(pdir, sprintf('%s_lev%d.dat', method_tags{im}, max_levels(il))),'w');
    fprintf(fid_p,'r\tst_h\tst_e\tsr_h\n');
    fprintf(fid_p,'%.6f\t%.8e\t%.8e\t%.8e\n', Pm.');
    fclose(fid_p);
    M.ndof=chsp{end}.ndof; M.nel=chm{end}.nel; M.wall=wall;
    R.(key)=M;
    fprintf('  ndof=%d nel=%d | Eg=%.3e Eb=%.3e Ei=%.2f%% Es=%.3e J=%.3f V=%.3e TV=%.3f (%.0fs)\n',...
        M.ndof,M.nel,M.E_glob,M.E_band,M.E_inf,M.E_slope,M.J,M.V_yield,M.TV,wall);
    save(fullfile(results_dir,'metrics.mat'),'R','method_tags','max_levels','c','dband','djump');
end
end

write_tables(R, method_tags, max_levels, results_dir, 'sphere');
fprintf('\nAll done. Results in %s\n', results_dir);

%% ---- exact Hill sphere hoop stress and its radial derivative ----
function [st,dst] = hill_sph_t(r,b,c,sy)
    st=zeros(size(r)); dst=zeros(size(r));
    ip = r <= c;
    st(ip)  = 2*sy*(0.5 - log(c./r(ip)) - (1-c^3/b^3)/3);
    dst(ip) = 2*sy./r(ip);
    st(~ip)  = 2*sy*c^3/(3*b^3)*(0.5*b^3./r(~ip).^3 + 1);
    dst(~ip) = -sy*c^3./r(~ip).^4;
end
