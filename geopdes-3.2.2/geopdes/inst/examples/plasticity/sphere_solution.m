% close all
% clear all
% clc

function [u, pressure] = sphere_solution(E,v, s_y, a, b)
% E =210000;
% v = 0.3;
% s_y = 0.24;
% a = 100;
% b = 200;

P_lim = 2 *s_y * log(b/a);
P_0 = 2*s_y/3*(1-(a^3/b^3));

P = linspace(0, P_lim*.9, 100);
ub = zeros(size(P));

% initial branch
for i =1:length(P)

        if P(i) < P_0
            ub(i) = 3*P(i) *b / ( 2*E* (b^3/a^3-1))*(1-v);

        else
            fun_front = @(x) -P(i) +2*s_y* log(x/a) + 2/3*s_y* (1- x.^3/b^3);
            c = fsolve(fun_front, 100);
            ub(i) = s_y*c^3/(E*b^2) * (1-v);

        end
end

% draw plateau
if P(end) >= P_0
    uc = linspace(ub(end), ub(end)*2, 100);
    F = zeros(size(uc));
    
    for i =1:length(uc)
        
        c = (uc(i) * E*b^2/(s_y* (1-v)))^(1/3);
        if c > b
            c = b;
        end
        F(i) = 2*s_y * log(c/a) + 2/3*s_y *(1- c^3/b^3);
    
    end
    
    u = [ub(1:end-1), uc];
    pressure = [P(1:end-1), F];
else
    u = ub;
    pressure = P;
end

% 
% figure(100)
% plot(u, pressure)

end
