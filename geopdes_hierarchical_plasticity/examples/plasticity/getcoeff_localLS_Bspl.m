% GETCOEFF_LOCALLS_BSPL: computes hierarchical quasi-interpolation coefficients 
% by (penalized) least squares in the local tensor-product spline space
%
% INPUT    
% hspace:    hierarchical space
% hmsh:    hierarchical mesh
% data: Mxd1 matrix containing coordinates of the points in the parametric space
% f: Mxd2 matrix containing values of the function to be approximated at the points in data
% lambda:  coefficient of the penalization term
% weight: weights of the points
%
% OUTPUT
% QI_coeff: coefficients of the solution
% 
% Author: Cesare Bracco

function QI_coeff = getcoeff_localLS_Bspl(hspace,hmsh,data,f,lambda,weight)

if nargin<6
    [M,~] = size(data);
    weight=ones(1,M);
end

par_dim=length(hspace.space_of_level(1).degree); %Dimension of parametric domain
deg=hspace.space_of_level(1).degree; %vector of degrees of the spline space
maxdeg=min([hspace.space_of_level(1).degree]); %maximum degree of the polynomial approximation 

%Auxiliary functions to handle cell-arry variables
extract_range=@(x,range) x(range);
last=@(x) max(x)+1;


lev=1;
if hspace.ndof_per_level(lev)~=0
    mass_matrix_full=op_gradgradu_gradgradv_tp(hspace.space_of_level(lev), hspace.space_of_level(lev), hmsh.mesh_of_level(lev));
    col_matrix_full=basisfun_multi(deg, data, hspace.space_of_level(lev).knots);
end
ndof_prev_levs=0;

for k=1:hspace.ndof
    %assemblaggio matrice
    while hspace.ndof_per_level(lev)==0 && lev<hspace.nlevels
        lev=lev+1;
        if hspace.ndof_per_level(lev)~=0
            mass_matrix_full=op_gradgradu_gradgradv_tp(hspace.space_of_level(lev), hspace.space_of_level(lev), hmsh.mesh_of_level(lev));
            col_matrix_full=basisfun_multi(deg, data, hspace.space_of_level(lev).knots);
        end     
    end
    kl=hspace.active{lev}(k-ndof_prev_levs); %index of current function in the B-spline basis of its level
     
    k_vec=cell(1,par_dim);
    [k_vec{:}]=ind2sub ([hspace.space_of_level(lev).ndof_dir, 1], kl); 
    support=sp_get_cells(hspace.space_of_level(lev),hmsh.mesh_of_level(lev),kl); %indices of the cells of the support of the B-spline
    I_vec=cell(1,par_dim);
    [I_vec{:}]=ind2sub([hmsh.mesh_of_level(lev).nel_dir, 1], support); 
    %Indices of the knots which are the corners of the support (the indices you get are numbered starting from 1)
    mu_vec=cellfun(@min,I_vec);
    nu_vec=cellfun(last,I_vec);
    
    %Corresponding coordinates
    breaks=hmsh.mesh_of_level(lev).breaks;
    p_mu=cellfun(extract_range,breaks,num2cell(mu_vec));
    p_nu=cellfun(extract_range,breaks,num2cell(nu_vec));
    
    nel_supp=nu_vec-mu_vec;
    nloc_min=prod(nel_supp+deg); %minimum local number of points that we require
    
    %Local set of data points (at first inside support)
    data_loc_ind=find(data(:,1)>=p_mu(1) & data(:,1)<=p_nu(1));
    for i=2:par_dim
        data_loc_ind=intersect(data_loc_ind,find(data(:,i)>=p_mu(i) & data(:,i)<=p_nu(i)));
    end 
    
    while numel(data_loc_ind)<nloc_min && numel(data_loc_ind)<length(f)
        mu_vec=max([mu_vec-1;ones(1,length(mu_vec))]);
        nu_vec=min([nu_vec+1;hmsh.mesh_of_level(lev).nel_dir+1]);
        p_mu=cellfun(extract_range,breaks,num2cell(mu_vec));
        p_nu=cellfun(extract_range,breaks,num2cell(nu_vec));
        data_loc_ind=find(data(:,1)>=p_mu(1) & data(:,1)<=p_nu(1));
        for i=2:par_dim
            data_loc_ind=intersect(data_loc_ind,find(data(:,i)>=p_mu(i) & data(:,i)<=p_nu(i)));
        end 
        support=max([support(1)-prod(hmsh.mesh_of_level(lev).nel_dir(1:end-1))-1 1]):...
        min([support(end)+prod(hmsh.mesh_of_level(lev).nel_dir(1:end-1))+1 prod(hmsh.mesh_of_level(lev).nel_dir(1:end))]);
    end 
    %Indices of B-splines active on local set of points
    ind_active_on_supp=sp_get_basis_functions(hspace.space_of_level(lev),...
    hmsh.mesh_of_level(lev), support);

    ind_loc_point{k}=data_loc_ind';  

   if numel(ind_loc_point{k})<(maxdeg+1)*(maxdeg+2)/2
        numel(ind_loc_point{k})
        error('Locally not enough data')
   end
    
    %Solving local (penalized) least squares approximation problem (with B-splines)
    loc_weight=weight(ind_loc_point{k});
    if max(ind_active_on_supp)>size(mass_matrix_full)
        keyboard
    end
    mass_matrix=mass_matrix_full(ind_active_on_supp,ind_active_on_supp);
    %col_matrix=basisfun_multi(deg, data(ind_loc_point{k},:), hspace.space_of_level(lev).knots); %collocation matrix for B-spline basis 
    col_matrix=col_matrix_full(ind_active_on_supp,ind_loc_point{k});
    bvector=f(ind_loc_point{k},:);
    bbbvector=col_matrix*diag(loc_weight)*bvector(:,:);
    S=col_matrix*diag(loc_weight)*col_matrix'+ lambda*mass_matrix; %lambda scalar
    loc_ind_Bspl_k=find(ind_active_on_supp==kl);
    for ivar = 1:size(f,2)
        QI_coeff_loc=S\bbbvector(:,ivar);         
        QI_coeff(k,ivar)=QI_coeff_loc(loc_ind_Bspl_k);
    end
    
    while k>=ndof_prev_levs+hspace.ndof_per_level(lev) && lev<hspace.nlevels
        ndof_prev_levs=ndof_prev_levs+hspace.ndof_per_level(lev);
        lev=lev+1;
        if hspace.ndof_per_level(lev)~=0
            mass_matrix_full=op_gradgradu_gradgradv_tp(hspace.space_of_level(lev), hspace.space_of_level(lev), hmsh.mesh_of_level(lev));
            col_matrix_full=basisfun_multi(deg, data, hspace.space_of_level(lev).knots);
        end 
    end
   
end
end