% close all
% clear all
% clc

function [u, pressure, sigma_r, sigma_t, radius] = Hill_solution(E,v, s_y, a, b, P_max, nload)
% E =210000;
% v = 0.3;
% s_y = 0.24;
% a = 100;
% b = 200;

P_lim = 2 *s_y/(sqrt(3)) * log(b/a);

Y = 2*s_y/(sqrt(3));
P_0 = Y/2*(1-(a^2/b^2));

P = linspace(P_max/nload, P_max, nload);

radius = linspace(a,b,1000);
ub = zeros(size(P));
sigma_r =  zeros(size(P,2), size(radius,2));
sigma_t =  zeros(size(P,2), size(radius,2));

% initial branch
for i =1:length(P)

        if P(i) < P_0
            ub(i) = 2*P(i)*b/(E*(b^2/a^2 -1))* (1-v^2);
            for j =1:length(radius)
                sigma_r(i,j) = (P(i)*a^2)/(b^2-a^2) - a^2*b^2*P(i)/((b^2-a^2)*radius(j)^2);
                sigma_t(i,j) = (P(i)*a^2)/(b^2-a^2) + a^2*b^2*P(i)/((b^2-a^2)*radius(j)^2);
            end
        else
            fun_front = @(x) -P(i)/Y + log(x/a) + .5* (1- x.^2/b^2);
            c = fsolve(fun_front, 100);
            ub(i) = Y*c^2/(E*b) * (1-v^2);

            for j =1:length(radius)

                if radius(j)<=c
                    sigma_r(i,j) = Y*(-.5 - log(c/radius(j) ) + c^2/(2*b^2));
                    sigma_t(i,j) = Y*(.5 - log(c/radius(j) ) + c^2/(2*b^2));
                else
                    sigma_r(i,j) = -Y * c^2/(2*b^2)*(b^2./(radius(j).^2) -1);
                    sigma_t(i,j) = Y * c^2/(2*b^2)*(b^2./(radius(j).^2) +1);
                end
            end
        end
end

% draw plateau
if P(end) >= P_0
    uc = linspace(ub(end), ub(end)*5, 100);
    F = zeros(size(uc));
    
    for i =1:length(uc)
        
        c = sqrt(uc(i) * E*b/(Y* (1-v^2)));
        if c > b
            c = b;
        end
        F(i) = Y * (log(c/a) + .5* (1- c^2/b^2));
    
    end
    
    u = [ub(1:end-1), uc];
    pressure = [P(1:end-1), F];
else
    u = ub;
    pressure = P;
end



end
