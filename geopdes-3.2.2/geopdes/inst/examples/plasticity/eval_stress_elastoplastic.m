function [internal_energy, K_tg, eps_pl, sigma_tot] = eval_stress_elastoplastic(eps_tot, eps_pl, mu, kappa, sigma_y, epsu, epsv, jacw, ncomp)
%EVAL_STRESS_ELASTOPLASTIC 
% perfect plasticity Neto, pages 221-222

if ncomp == 2
    eps_tot = [eps_tot(1), eps_tot(4), 0, eps_tot(2), 0, 0]';
elseif ncomp == 3
    eps_tot = [eps_tot(1), eps_tot(5), eps_tot(9), eps_tot(2),  eps_tot(3),  eps_tot(6)]';
end


[sigma_tot, C_tg, eps_pl] = perfect_plasticity_model(eps_tot, eps_pl, mu, kappa, sigma_y);


if ncomp == 2
    epsv = reshape(epsv, [size(epsv,1),size(epsv,3)]);
    epsu = reshape(epsu, [size(epsu,1),size(epsu,4)]);
    epsv = [epsv(1,:); epsv(4,:); 2*epsv(2,:)];
    epsu = [epsu(1,:); epsu(4,:); 2*epsu(2,:)];

    sigma = [sigma_tot(1) sigma_tot(2) sigma_tot(4)]';
    C_tg = [C_tg(1,1) C_tg(1,2) C_tg(1,4); C_tg(2,1) C_tg(2,2) C_tg(2,4); C_tg(4,1) C_tg(4,2) C_tg(4,4) ];
elseif ncomp == 3
    sigma = sigma_tot;
    epsv = reshape(epsv, [size(epsv,1),size(epsv,3)]);
    epsu = reshape(epsu, [size(epsu,1),size(epsu,4)]);
    epsv = [epsv(1,:); epsv(5,:); epsv(9,:); 2*epsv(2,:); 2*epsv(3,:); 2*epsv(6,:)];
    epsu = [epsu(1,:); epsu(5,:); epsu(9,:); 2*epsu(2,:); 2*epsu(3,:); 2*epsu(6,:)];
end

internal_energy = epsv' * sigma * jacw;                     % internal energy
K_tg = epsv' * C_tg * epsu * jacw;                          % tangent matrix

end

