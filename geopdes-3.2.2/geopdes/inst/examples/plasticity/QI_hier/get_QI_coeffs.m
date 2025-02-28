function [ QI_coeff, ind_loc_point, muxs,nuxs,muys,nuys, condition, all_loc,local_app,min_max_msv_basischng] = get_QI_coeffs(hspace,hmsh,data)
% INPUT:
% hspace: hierarchical space
% hmsh: hierarchcial space
% data: n x(2+m) matrix, where n is the number of data and m is the vector dimension of the data values
% OUTPUT: 
% QI_coeff: (hspace.ndof) x m
% ind_loc_points: cell-array containing, for each basis function, the
% indices of the data points used to compute its coefficient
% other technical outputs that we do not need now...

f=data(:,3:end);
data=data(:,1:2);
a=min(data(:,1));
b=max(data(:,1));
c=min(data(:,2));
d=max(data(:,2));
data(:,1)=(data(:,1)-a)/(b-a);
data(:,2)=(data(:,2)-c)/(d-c);
%cell array of hspace.ndof_per_level(i)x2
local_app=cell(hspace.nlevels,1);
polynomial=0; % flag to switch to polynomial local least squares
dx=hspace.space_of_level(1).degree(1);
dy=hspace.space_of_level(1).degree(2);
soglia=10^2; %blackforest 3*10^3 %MTU tests 10^8 %SOGLIA
maxdeg=dx;

element=cell(hspace.ndof,4);
lev=1;
ndof_prev_levs=0;
min_max_msv_basischng(lev,1)=10^20;
min_max_msv_basischng(lev,2)=0;

for k=1:hspace.ndof
    if (k==1)
        kprec=0;
    else
        kprec=length(QI_coeff);
    end
    while hspace.ndof_per_level(lev)==0 && lev<hspace.nlevels
        lev=lev+1;
        min_max_msv_basischng(lev,1)=10^20;
        min_max_msv_basischng(lev,2)=0;
    end
    kl=hspace.active{lev}(k-ndof_prev_levs);
    curr_nx=hmsh.mesh_of_level(lev).nel_dir(1);
    curr_ny=hmsh.mesh_of_level(lev).nel_dir(2);
    %clamped knot vectors
    x_clamp=augknt(linspace(0,1,curr_nx+1),dx+1);
    y_clamp=augknt(linspace(0,1,curr_ny+1),dy+1);
    
    [k_x,k_y]=ind2sub (hspace.space_of_level(lev).ndof_dir, kl);    
    support=sp_get_cells(hspace.space_of_level(lev),hmsh.mesh_of_level(lev),(k_y-1)*(curr_nx+dx)+k_x); %indices of the cells of the support of the B-spline
    [I,J]=ind2sub(hmsh.mesh_of_level(lev).nel_dir, support); 
    %Indices of the knots which are the corners of the support (the indices you get are numbered starting from 1)
    I1_x=min(I);
    I2_x=max(I)+1;
    I1_y=min(J);
    I2_y=max(J)+1;
    %Indices \mu and \nu delimiting the local subset which determines the B-splines involved in the local linear system
    mux=I1_x;
    muy=I1_y;

    nux=I2_x;
    nuy=I2_y;
    element{k}=[mux nux muy nuy];
    
    %Corresponding abscissae and ordinates
    breaksx=hmsh.mesh_of_level(lev).breaks{1};
    breaksy=hmsh.mesh_of_level(lev).breaks{2};
    x_mux=0+breaksx(mux);%0+(mux-1)/curr_nx; uniform case
    x_nux=0+breaksx(nux);%0+(nux-1)/curr_nx;
    y_muy=0+breaksy(muy);%0+(muy-1)/curr_ny;
    y_nuy=0+breaksy(nuy);%0+(nuy-1)/curr_ny;
    
    %First we attempt to get the coefficient by a local least squares
    %approximation in the spline space
    data_loc_x=find(data(:,1)>=x_mux & data(:,1)<=x_nux);
    data_loc_y=find(data(:,2)>=y_muy & data(:,2)<=y_nuy);
    data_loc_ind=intersect(data_loc_x,data_loc_y); 
    all_loc{k}=data_loc_ind;  
    ind_loc_point{k}=data_loc_ind';
    
    polynomial=1;
    %radius of the circle containing the support of the B-spline
    %Beginning of variant
%     if (mux-nux)==dx+1
%         x_mux=x_mux+1/(2*curr_nx);
%         x_nux=x_nux-1/(2*curr_nx);
%     end
%     if (muy-nuy)==dy+1
%         y_muy=y_muy+1/(2*curr_ny);
%         y_nuy=y_nuy-1/(2*curr_ny);
%     end
    %end of variant
    r=sqrt((x_nux-x_mux)^2+(y_nuy-y_muy)^2)/2;
    if lev>=2
        ratioray=ceil(sqrt(((dx+1)/hmsh.mesh_of_level(lev-1).nel_dir(1))^2+((dy+1)/hmsh.mesh_of_level(lev-1).nel_dir(2))^2)/(2*r));
    else
        ratioray=ceil(sqrt(((dx+1)/hmsh.mesh_of_level(lev).nel_dir(1))^2+((dy+1)/hmsh.mesh_of_level(lev).nel_dir(2))^2)/(2*r));    
    end
    %vector with the coordinates of the center repeated
    center=ones(size(data,1),1)*[(x_mux+x_nux)/2 (y_muy+y_nuy)/2];
    data_dist=sqrt(sum((abs(data-center)).^2,2));
    data_loc_ind=find(data_dist<=r); 
    all_loc{k}=data_loc_ind;  
    ind_loc_point{k}=data_loc_ind';
    
    pd=maxdeg;
    rsaved=r;

    i=1;
    while (numel(ind_loc_point{k})<(pd+1)*(pd+2)/2) && (i<=2*ratioray+1) 
        r=i*rsaved;
        data_loc_ind=find(data_dist<=r); 
        all_loc{k}=data_loc_ind;  
        ind_loc_point{k}=data_loc_ind';
        i=i+1;
        pd=maxdeg;
        while (numel(ind_loc_point{k})<(pd+1)*(pd+2)/2) && (pd>0)
            pd=pd-1;
        end
    end            

    
   if numel(ind_loc_point{k})<(pd+1)*(pd+2)/2
        pd
        i
        numel(ind_loc_point{k})
        error('Locally not enough data')
   end
    
    index=0;
    fullrank=0;
    max_pd=pd;
    clear pd;
    for pd=max_pd:-1:0
        
    %solution of the polynomial least sqaures problem
    monx=eye(pd+1);
    mony=eye(pd+1);
    deg_vectorx=repelem([0:pd],pd+1);
    deg_vectory=repmat([0:pd],1,pd+1);
    sum_deg=deg_vectorx+deg_vectory;
    deg_vector=find(sum_deg<=pd);
    dataloc = data(ind_loc_point{k},:)';
    Dx = x_nux-x_mux;
    Dy = y_nuy-y_muy;
    %size(dataloc)
    datascal(1,:) = (dataloc(1,:)-x_mux)/Dx;
    datascal(2,:) = (dataloc(2,:)-y_muy)/Dy;
    for h=1:length(deg_vector)
        Blocpol_lsmatrix(h,:)=polyval(monx(pd+1-deg_vectorx(deg_vector(h)),:),datascal(1,:)).*polyval(mony(pd+1-deg_vectory(deg_vector(h)),:),datascal(2,:));
    end

    coeffs_polyloc=pinv(Blocpol_lsmatrix')*f(ind_loc_point{k},:);
    % riordino dei coefficienti per facilitare cambiamento di base
       switch pd
           case 0 % degree 0
                coeffs_poly = coeffs_polyloc;
                C=1; D=1;
           case 1 % degree 1
                coeffs_polyloc = coeffs_polyloc([1 3 2],:);
                
                D=diag([1 1/Dx,1/Dy]);
                x0 = x_mux; y0= y_muy;
                C= eye(3);
                C(2,1)=-x0; 
                C(3,1)= -y0;
                
                coeffs_poly = C'*D*coeffs_polyloc;
               
                coeffs_poly = coeffs_poly([1 3 2],:);
           case 2 % degree 2
                coeffs_polyloc = coeffs_polyloc([1 4 2 6 5 3],:); 
        
                D=diag([1 1/Dx,1/Dy,1/Dx^2,1/(Dx*Dy),1/Dy^2]);
                x0 = x_mux; y0= y_muy;
                C= eye(6);
                C(2,1)=-x0; 
                C(3,1)= -y0;
                C(4,1) = x0^2; C(4,2)=-2*x0;
                C(5,1)= x0*y0; C(5,2)=-y0; C(5,3) = -x0;
                C(6,1)= y0^2;  C(6,3) = -2*y0;
                coeffs_poly = C'*D*coeffs_polyloc;
                
                coeffs_poly = coeffs_poly([1 3 6 2 5 4],:);
           case 3
               coeffs_polyloc = coeffs_polyloc([1 5 2 8 6 3 10 9 7 4],:); 
        
                D=diag([1 1/Dx,1/Dy,1/Dx^2,1/(Dx*Dy),1/Dy^2,1/Dx^3,1/(Dx^2*Dy),1/(Dx*Dy^2),1/Dy^3]);
                x0 = x_mux; y0= y_muy;
                C= eye(10);
                C(2,1)=-x0; 
                C(3,1)= -y0;
                C(4,1) = x0^2; C(4,2)=-2*x0;
                C(5,1)= x0*y0; C(5,2)=-y0; C(5,3) = -x0;
                C(6,1)= y0^2;  C(6,3) = -2*y0;
                C(7,1)=-x0^3; C(7,2)=3*x0^2; C(7,4)=-3*x0;
                C(8,1)=-x0^2*y0; C(8,2)=2*x0*y0; C(8,3)=x0^2; C(8,4)=-y0; C(8,5)=-2*x0;
                C(9,1)=-x0*y0^2; C(9,2)=y0^2; C(9,3)=2*x0*y0; C(9,5)=-2*y0; C(9,6)=-x0;
                C(10,1)=-y0^3; C(10,3)=3*y0^2; C(10,6)=-3*y0;
                coeffs_poly = C'*D*coeffs_polyloc;
                
                coeffs_poly = coeffs_poly([1 3 6 10 2 5 9 4 8 7],:);
           case 4
               coeffs_polyloc = coeffs_polyloc([1 6 2 10 7 3 13 11 8 4 15 14 12 9 5],:); 
        
                D=diag([1 1/Dx,1/Dy,1/Dx^2,1/(Dx*Dy),1/Dy^2,1/Dx^3,1/(Dx^2*Dy),1/(Dx*Dy^2),1/Dy^3,...
                    1/Dx^4 1/(Dx^3*Dy) 1/(Dx^2*Dy^2) 1/(Dx*Dy^3) 1/(Dy^4)]);
                x0 = x_mux; y0= y_muy;
                C= eye(15);
                C(2,1)=-x0; 
                C(3,1)= -y0;
                C(4,1) = x0^2; C(4,2)=-2*x0;
                C(5,1)= x0*y0; C(5,2)=-y0; C(5,3) = -x0;
                C(6,1)= y0^2;  C(6,3) = -2*y0;
                C(7,1)=-x0^3; C(7,2)=3*x0^2; C(7,4)=-3*x0;
                C(8,1)=-x0^2*y0; C(8,2)=2*x0*y0; C(8,3)=x0^2; C(8,4)=-y0; C(8,5)=-2*x0;
                C(9,1)=-x0*y0^2; C(9,2)=y0^2; C(9,3)=2*x0*y0; C(9,5)=-2*y0; C(9,6)=-x0;
                C(10,1)=-y0^3; C(10,3)=3*y0^2; C(10,6)=-3*y0;
                C(11,1)=x0^4; C(11,2)=-4*x0^3; C(11,4)=6*x0^2; C(11,7)=-4*x0;
                C(12,1)=x0^3*y0; C(12,2)=-3*x0^2*y0; C(12,3)=-x0^3; C(12,4)=3*x0*y0; C(12,5)=3*x0^2; C(12,7)=-y0; C(12,8)=-3*x0;
                C(13,1)=x0^2*y0^2; C(13,2)=-2*x0*y0^2; C(13,3)=-2*x0^2*y0; C(13,4)=y0^2; C(13,5)=4*x0*y0; C(13,6)=x0^2; C(13,8)=-2*y0; C(13,9)=-2*x0;
                C(14,1)=x0*y0^3; C(14,2)=-y0^3; C(14,3)=-3*x0*y0^2; C(14,5)=3*y0^2; C(14,6)=3*x0*y0; C(14,10)=-x0; C(14,9)=-3*y0;
                C(15,1)=y0^4; C(15,3)=-4*y0^3; C(15,6)=6*y0^2; C(15,10)=-4*y0;
                coeffs_poly = C'*D*coeffs_polyloc;
                
                coeffs_poly = coeffs_poly([1 3 6 10 15 2 5 9 14 4 8 13 7 12 11],:);
       end
        %computation of the coefficients of the polynomial in the B-basis
    x_est=hspace.space_of_level(lev).knots{1}; %uniform case x_est=augknt(linspace(0,1,curr_nx+1),dx+1);%uniform case
    y_est=hspace.space_of_level(lev).knots{2}; %y_est=augknt(linspace(0,1,curr_ny+1),dy+1);
    col_pointsx=linspace(x_mux,x_nux,dx+nux-mux);
    col_pointsy=linspace(y_muy,y_nuy,dx+nuy-muy);
    aux_x=(spcol(x_est,dx+1,col_pointsx))';
    aux_y=(spcol(y_est,dy+1,col_pointsy))';
    for h=mux:nux+dx-1
        for q=muy:nuy+dy-1    
            aux=aux_x(h,:)'*aux_y(q,:);
            bbasis_lmatrix((h-mux)*(dy+nuy-muy)+q-muy+1,:)=aux(:)';
        end
    end
    %auxiliary matrix to evaluate on the grid the bivariate polynomial
    %solution of the least squares problem
    for h=1:length(deg_vector)  
        auxp=(polyval(monx(pd+1-deg_vectorx(deg_vector(h)),:),col_pointsx))'*(polyval(mony(pd+1-deg_vectory(deg_vector(h)),:),col_pointsy));
        poly_ematrix(h,:)=auxp(:)';
    end 
    b=poly_ematrix'*coeffs_poly;
    coeffs_BB=(bbasis_lmatrix')\b;
    
    A_J_pseudoinverse=pinv(Blocpol_lsmatrix');  
    G_J=inv(bbasis_lmatrix')*poly_ematrix'*C'*D*pinv(Blocpol_lsmatrix');
    T1 = svd(G_J);
    MSV = max(abs(T1));
    T1_2 = svd(A_J_pseudoinverse);
    MSV2= max(abs(T1_2));
    if MSV2<soglia
    index=index+1;
    adm_deg(index)=pd;
    adm_MSV(index)=MSV2;
    
    ind_coeffx=find([mux:nux+dx-1]==k_x);
    ind_coeffy=find([muy:nuy+dy-1]==k_y);
    ind_coeff=(ind_coeffx-1)*(dy+nuy-muy)+ind_coeffy;
    adm_QI_coeff(index,:)=coeffs_BB(ind_coeff,:);
    
    %saving minimum and maximum value of maximal singular vakues of the
    %basis change matrix
    min_max_msv_basischng(lev,1)=min([min_max_msv_basischng(lev,1),max(abs(svd(inv(bbasis_lmatrix')*poly_ematrix')))]);
    min_max_msv_basischng(lev,2)=max([min_max_msv_basischng(lev,2),max(abs(svd(inv(bbasis_lmatrix')*poly_ematrix')))]);
    
    end
    
    clear Blocpol_lsmatrix
    clear coeffs_polyloc
    clear datascal
    %clear polynomial_lsmatrix
    clear coeffs_poly
    clear poly_ematrix
    clear b
    clear bbasis_lmatrix
    clear coeffs_BB
    end

   [pd, choice]=max(adm_deg);
    QI_coeff(k,:)=adm_QI_coeff(choice,:);
    
    clear adm_MSV
    clear adm_deg
    clear adm_QI_coeff
    
    
    local_app{lev}(kl,1)=polynomial;
    if polynomial==1
        local_app{lev}(kl,2)=max(0,pd);
    end
    while k>=ndof_prev_levs+hspace.ndof_per_level(lev) && lev<hspace.nlevels
        ndof_prev_levs=ndof_prev_levs+hspace.ndof_per_level(lev);
        lev=lev+1;
        min_max_msv_basischng(lev,1)=10^20;
        min_max_msv_basischng(lev,2)=0;
    end
    muxs(k)=mux;
    nuxs(k)=nux;
    muys(k)=muy;
    nuys(k)=nuy;
    condition=0;
end
end