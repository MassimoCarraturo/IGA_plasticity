function est = adaptivity_estimate_distance_plastic_front(hmsh,problem_data,mult)

nel = hmsh.nel;
est = zeros(nel,1);

% position plastic front

a= problem_data.R_i; % inner radius
b = problem_data.R_o; % outer radius
s_y = problem_data.s_y; %yielding stress
P_max = problem_data.Pmax;
P = P_max*mult;


Y = 2*s_y/(sqrt(3));
P_0 = Y/2*(1-(a^2/b^2));
if P<= P_0
    R = a;
else
    fun_front = @(x) -P/Y + log(x/a) + .5* (1- x.^2/b^2);
    R = fsolve(fun_front, 100);     
end

% compute distance for each element
k = 1;
for ilev =1:hmsh.nlevels %loop levels
    nel_act_lev = numel(hmsh.active{ilev});
    if nel_act_lev >0
        for iel =1:nel_act_lev % loop active elements of level
            % id_el = hmsh.active{ilev}(iel);
            coo_gp = hmsh.msh_lev{ilev}.geo_map(:,:,iel); % gauss point coordinates
            nqn = size(coo_gp,2);
            avg_dist = 0; % average distance of the gauss points from the plastic front
            for ipt =1:nqn
                r = sqrt(coo_gp(1,ipt)^2 + coo_gp(2,ipt)^2);
                avg_dist = avg_dist + abs(r -R);
            end
            avg_dist = avg_dist/nqn;
            est(k) = - avg_dist;
            k = k+1;
        end
    end
end

min_est = min(est);
est = est - min_est;

end