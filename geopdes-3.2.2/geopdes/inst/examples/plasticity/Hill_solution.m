% close all
% clear all
% clc

function [u, pressure] = Hill_solution(E,v, s_y, a, b)
% E =210000;
% v = 0.3;
% s_y = 0.24;
% a = 100;
% b = 200;

P_lim = 2 *s_y/(sqrt(3)) * log(b/a);

Y = 2*s_y/(sqrt(3));
P_0 = Y/2*(1-(a^2/b^2));

P = linspace(0, P_lim*.9, 100);
ub = zeros(size(P));

for i =1:length(P)

        if P(i) < P_0
            ub(i) = 2*P(i)*b/(E*(b^2/a^2 -1))* (1-v^2);
        else
            fun_front = @(x) -P(i)/Y + log(x/a) + .5* (1- x.^2/b^2);
            c = fsolve(fun_front, 100);
            ub(i) = Y*c^2/(E*b) * (1-v^2);

        end


end

uc = linspace(ub(end), .6, 100);
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
% figure(100)
% plot(u, pressure)

end
