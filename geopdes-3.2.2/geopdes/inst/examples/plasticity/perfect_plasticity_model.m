function [sigma, C_tg, eps_pl] = perfect_plasticity_model(eps_tot, eps_pl, mu, kappa, sigma_y)


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

end