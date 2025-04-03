function est =  adaptivity_estimate_div_sigma_el (sigma_store, hmsh, hspace, hspace_scalar, problem_data, adaptivity_data);

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

divergence = zeros(hspace.ncomp, hmsh.mesh_of_level.nqn, hmsh.nel);
for ifield =1:6
    [hgrad, F] = hspace_eval_hmsh (sigma_store(:,ifield), hspace_scalar, hmsh, 'gradient');

    if ifield ==1
        divergence(1,:,:) = divergence(1,:,:) + hgrad(1,:,:); %sigma_xx,x
    elseif ifield ==2
        divergence(2,:,:) = divergence(2,:,:) + hgrad(2,:,:); %sigma_yy,y           
    elseif ifield ==4
        divergence(1,:,:) = divergence(1,:,:) + hgrad(2,:,:); %sigma_xy,y
        divergence(2,:,:) = divergence(2,:,:) + hgrad(1,:,:); %sigma_yx,x

    if hspace.ncomp ==3
        elseif ifield ==3
            divergence(3,:,:) = divergence(3,:,:) + hgrad(3,:,:); %sigma_zz,z
        elseif ifield ==5
            divergence(1,:,:) = divergence(1,:,:) + hgrad(3,:,:) ;%sigma_xz,z
            divergence(3,:,:) = divergence(3,:,:) + hgrad(1,:,:); %sigma_zx,x
        elseif ifield ==6    
            divergence(2,:,:) = divergence(2,:,:) + hgrad(3,:,:); %sigma_yz,z
            divergence(3,:,:) = divergence(3,:,:) + hgrad(2,:,:); %sigma_zy,y
        end
    end 

end



x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim  %rdim is the dimension of the physical domain
    x{idim} = reshape (F(idim,:), [], hmsh.nel);
end




valf = problem_data.f (x{:}) ;
% for h = 1:hspace.ncomp
%     divergence(h,:,:) = 
% end
aux = (valf + divergence).^2;  %residual

switch adaptivity_data.flag
    case 'elements'
        w = [];
        h = [];
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
                h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
            end
        end
        h = h * sqrt (hmsh.ndim);

        aux = reshape (sum (aux), [], hmsh.nel);
        est = sum (aux.*w);
        est = h.^2 .* est(:);

        % Jump terms, only computed for multipatch geometries
        % NOT YET IMPLEMENTED
        est = C0_est * sqrt (est);

    case 'functions'
        ms = zeros (hmsh.nlevels, 1);
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
            else
                ms(ilev) = 0;
            end
        end
        ms = ms * sqrt (hmsh.ndim);
        
        Nf = cumsum ([0; hspace.ndof_per_level(:)]);
        dof_level = zeros (hspace.ndof, 1);
        for lev = 1:hspace.nlevels
            dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
        end
        coef = ms(dof_level).^2 .* hspace.coeff_pou(:);
        
        est = zeros(hspace.ndof,1);
        ndofs = 0;
        Ne = cumsum([0; hmsh.nel_per_level(:)]);
        for ilev = 1:hmsh.nlevels
            ndofs = ndofs + hspace.ndof_per_level(ilev);
            if (hmsh.nel_per_level(ilev) > 0)
                ind_e = (Ne(ilev)+1):Ne(ilev+1);
                sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
                b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, aux(:,:,ind_e));
                dofs = 1:ndofs;
                est(dofs) = est(dofs) + hspace.Csub{ilev}.' * b_lev;
            end
        end
        est = coef .* est;

        % Jump terms, only computed for multipatch geometries 
        % NOT YET IMPLEMENTED
        est = C0_est * sqrt (est);
end

end