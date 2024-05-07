function [internal_energy, K_tg, eps_pl] = eval_stress_elastoplastic(eps_tot, eps_pl, mu, kappa, sigma_y, epsu, epsv, jacw, ncomp)

%EVAL_STRESS_ELASTOPLASTIC 
if ncomp == 2
    eps_tot = [eps_tot(1), eps_tot(4), 0, eps_tot(2), 0, 0]';
elseif ncomp == 3
    eps_tot = [eps_tot(1), eps_tot(5), eps_tot(9), eps_tot(2),  eps_tot(3),  eps_tot(6)]';
end

theta = [1 1 1 0 0 0]*eps_tot;

eps_tot_dev = eps_tot-1/3*[1 1 1 0 0 0]'*theta;         % deviatoric tital strain
sigma_dev_trial = 2*mu*(eps_tot_dev - eps_pl);          % deviatoric stress trial
sigm_dev_norm = sqrt(sigma_dev_trial(1)^2+sigma_dev_trial(2)^2+sigma_dev_trial(3)^2+...
    2*(sigma_dev_trial(4)^2+sigma_dev_trial(5)^2+sigma_dev_trial(6)^2));
q_trial = sqrt(3/2)*sigm_dev_norm;              % elastic trial von Mises effective stress

Idev = [2/3 -1/3 -1/3 0 0 0; -1/3 2/3 -1/3 0 0 0; -1/3 -1/3 2/3 0 0 0;...
    0 0 0 1/2 0 0; 0 0 0 0 1/2 0; 0 0 0 0 0 1/2];
C_tg_el = kappa*[ones(3), zeros(3); zeros(3,6)] + 2*mu*Idev;

if q_trial > sigma_y
    % evaluate sigma
    lamda = (q_trial-sigma_y)/(3*mu);                           % incremental plastic multiplier
    n_trial = sigma_dev_trial/sigm_dev_norm;                    % norm of the yield surface
    C = 3*mu*lamda/q_trial;
    sigma = sigma_dev_trial*(1-C)+ kappa*theta*[1 1 1 0 0 0]';  % stress tensor
        
    % evaluate tangent modulus 
    N = n_trial*n_trial';                                       % unit flow vector
    C_tg = C_tg_el -...                                         % elastoplastic consistent tangent modulus for von Mises model
        2*mu*C*Idev + ...
        2*mu*(C-1)*N;
    % update plastic strain
    eps_pl = eps_pl + sqrt(3/2)*n_trial*lamda;                  % update of the plastic strain tensor
else
    % evaluate sigma
    sigma= sigma_dev_trial + kappa*theta*[1 1 1 0 0 0]';
    % evaluate tangent
    C_tg = C_tg_el;
end

if ncomp == 2
    epsv = reshape(epsv, [size(epsv,1),size(epsv,3)]);
    epsu = reshape(epsu, [size(epsu,1),size(epsu,4)]);
    epsv = [epsv(1,:); epsv(4,:); 2*epsv(2,:)];
    epsu = [epsu(1,:); epsu(4,:); 2*epsu(2,:)];

    sigma = [sigma(1) sigma(2) sigma(4)]';
    C_tg = [C_tg(1,1) C_tg(1,2) C_tg(1,4); C_tg(2,1) C_tg(2,2) C_tg(2,4); C_tg(4,1) C_tg(4,2) C_tg(4,4) ];
elseif ncomp == 3
    epsv = reshape(epsv, [size(epsv,1),size(epsv,3)]);
    epsu = reshape(epsu, [size(epsu,1),size(epsu,4)]);
    epsv = [epsv(1,:); epsv(5,:); epsv(9,:); 2*epsv(2,:); 2*epsv(3,:); 2*epsv(6,:)];
    epsu = [epsu(1,:); epsu(5,:); epsu(9,:); 2*epsu(2,:); 2*epsu(3,:); 2*epsu(6,:)];
end

internal_energy = epsv' * sigma * jacw;                     % internal energy
K_tg = epsv' * C_tg * epsu * jacw;                          % tangent matrix

end

