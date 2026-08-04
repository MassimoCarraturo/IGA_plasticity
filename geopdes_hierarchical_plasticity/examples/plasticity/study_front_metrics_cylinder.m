% STUDY_FRONT_METRICS_CYLINDER  Front-focused metrics for the projection
% operators on the pressurised thick-walled cylinder (exact Hill reference).
%
% The mesh is driven by the GEOMETRIC front estimator, which is solution
% independent, so every operator runs on an IDENTICAL front-refined mesh and
% the metrics measure the projector alone.  Metric definitions: see
% front_kink_metrics.m.
%
% Exact kink: sigma_t is continuous at r = c but its radial derivative jumps,
%   plastic (r<c): sigma_t = Y(1/2 - ln(c/r) + c^2/(2 b^2)),  d/dr = +Y/r
%   elastic (r>c): sigma_t = Y c^2/(2 b^2)(b^2/r^2 + 1),      d/dr = -Y c^2/r^3
% so [d sigma_t/dr] = 2 Y / c.   (sigma_r is C^1 there.)
%
% Output: results/front_metrics_cylinder/

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
results_dir = fullfile(here, 'results', 'front_metrics_cylinder');
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

%% physical data
E=210000; nu=0.49; sigma_y=240; a=100; b=200; nload=5;
Y = 2*sigma_y/sqrt(3);
P = Y*( log(160/a) + 0.5*(1-160^2/b^2) );

pd.geo_name='geo_ring_SouzaNeto.txt';
pd.nmnn_sides=[2]; pd.drchlt_sides=[]; pd.press_sides=[1];
pd.symm_sides=[3 4]; pd.slider_sides=[];
pd.yield_stress=@(x,y) sigma_y*ones(size(x));
pd.kappa_lame=@(x,y) E/(3*(1-2*nu))*ones(size(x));
pd.mu_lame=@(x,y) E/(2*(1+nu))*ones(size(x));
pd.f=@(x,y) zeros(2,size(x,1),size(x,2));
pd.g=@(x,y,ind) zeros(2,size(x,1),size(x,2));
pd.h=@(x,y,ind,mult) zeros(2,size(x,1),size(x,2));
pd.p=@(x,y,ind) P*ones(size(x));
pd.R_i=a; pd.R_o=b; pd.s_y=sigma_y; pd.Pmax=P;

c = fzero(@(x) -P/Y + log(x/a) + .5*(1-x.^2/b^2), [a b]);
fprintf('front c = %.4f mm ; exact slope jump 2Y/c = %.4f MPa/mm\n', c, 2*Y/c);

%% discretization: geometric front marking -> identical mesh for all operators
p=4;
md0.degree=[p p]; md0.regularity=[p-1 p-1]; md0.nsub_coarse=[5 5]; md0.nsub_refine=[2 2];
md0.nquad=[p+1 p+1]; md0.space_type='standard'; md0.truncated=1; md0.nload=nload;
md0.newton_tol=1e-8; md0.newton_tol_abs=1e-10; md0.newton_iter_max=100;
ad0.flag='elements'; ad0.estimator='plastic_front'; ad0.C0_est=1.0;
ad0.mark_param=0.8; ad0.mark_param_coarsening=0.1; ad0.mark_strategy='MS';
ad0.max_ndof=1e7; ad0.max_nel=1e6; ad0.num_max_iter=3; ad0.tol=1e-10;   % must not bind: max_level is the control
ad0.adm_strategy='admissible'; ad0.coarsening_flag='any'; ad0.adm=p;

dband=15; djump=12; nsamp=6001;
method_tags={'L2','QI','QI_C0','QI_graded','QI_fine','Bezier','Bezier_C0','DLSQ','DLSQ_W'};
max_levels=[3 4];

R=struct();
for il=1:numel(max_levels)
for im=1:numel(method_tags)
    md=md0; ad=ad0; ad.max_level=max_levels(il);
    switch method_tags{im}
        case 'L2';        md.type_projection='L2';
        case 'QI';        md.type_projection='QI';
        case 'QI_C0';     md.type_projection='QI_C0';
        case 'QI_graded'; md.type_projection='QI_graded'; md.graded_transition_level=2;
        case 'QI_fine';   md.type_projection='QI'; md.num_bisections=1;
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
    uu=linspace(0,1,nsamp); vv=0.5;                 % theta = 45 deg ray
    [~,F]=sp_eval(sig(:,1),hs,geo,{uu,vv});
    r=sqrt(squeeze(F(1,:,:)).^2+squeeze(F(2,:,:)).^2); r=r(:);

    S=cell(1,6);
    for k=1:6, S{k}=reshape(sp_eval(sig(:,k),hs,geo,{uu,vv}),[],1); end
    st_h  = 0.5*(S{1}+S{2}) - S{4};                 % sigma_t on the 45 deg ray
    dst_h = 0.5*(gr(sig(:,1),hs,geo,uu,vv)+gr(sig(:,2),hs,geo,uu,vv)) - gr(sig(:,4),hs,geo,uu,vv);
    vm_h  = sqrt(0.5*((S{1}-S{2}).^2+(S{2}-S{3}).^2+(S{3}-S{1}).^2)+3*(S{4}.^2+S{5}.^2+S{6}.^2));

    [st_e,dst_e] = hill_cyl_t(r,b,c,Y);
    M = front_kink_metrics(r, st_h, dst_h, st_e, dst_e, vm_h, sigma_y, c, dband, djump);
    % Profile dump for the stress-kink figures: sigma_t and its analytic
    % counterpart along the sampling ray, subsampled to keep the files small.
    pdir = fullfile(results_dir,'profiles');
    if ~exist(pdir,'dir'); mkdir(pdir); end
    ss = max(1, round(numel(r)/800));
    idx = 1:ss:numel(r);
    Pm = [r(idx), st_h(idx), st_e(idx)];
    fid_p = fopen(fullfile(pdir, sprintf('%s_lev%d.dat', method_tags{im}, max_levels(il))),'w');
    fprintf(fid_p,'r\tst_h\tst_e\n');
    fprintf(fid_p,'%.6f\t%.8e\t%.8e\n', Pm.');
    fclose(fid_p);
    M.ndof=chsp{end}.ndof; M.nel=chm{end}.nel; M.wall=wall;
    R.(sprintf('%s_lev%d',method_tags{im},max_levels(il)))=M;
    fprintf('  ndof=%d nel=%d | Eg=%.3e Eb=%.3e Ei=%.2f%% Es=%.3e J=%.3f V=%.3e TV=%.3f\n',...
        M.ndof,M.nel,M.E_glob,M.E_band,M.E_inf,M.E_slope,M.J,M.V_yield,M.TV);
end
end

save(fullfile(results_dir,'metrics.mat'),'R','method_tags','max_levels','c','dband','djump');
write_tables(R, method_tags, max_levels, results_dir, 'cylinder');
fprintf('\nAll done. Results in %s\n', results_dir);

%% ---- helpers ----
function d = gr(coefs,hs,geo,uu,vv)
    g = sp_eval(coefs,hs,geo,{uu,vv},'gradient');
    d = (reshape(g(1,:,:),[],1)+reshape(g(2,:,:),[],1))/sqrt(2);
end

function [st,dst] = hill_cyl_t(r,b,c,Y)
    st=zeros(size(r)); dst=zeros(size(r));
    ip = r<c;
    st(ip)  = Y*(0.5-log(c./r(ip))+c^2/(2*b^2));   dst(ip)  = Y./r(ip);
    st(~ip) = Y*c^2/(2*b^2)*(b^2./r(~ip).^2+1);    dst(~ip) = -Y*c^2./r(~ip).^3;
end
