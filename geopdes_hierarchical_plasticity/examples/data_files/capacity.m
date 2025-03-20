function C = capacity(eu,eu_old)
%CAPACITY = volumetric heat capacity for Ti-6Al-4V from K. C. Mills,"Recommended Values of
%Thermophysic Properties for selected commercial alloys", 2002.
% 
% C = zeros(size(eu));
% 
% for i=1:numel(eu)
%     if eu(i) <= 950
%         density = 4420 - 0.154*(eu(i)-25);
%         density = density * 1e-9; % rescale to [kg/mm^3]
%         capacity = (0.55 + 0.2*(eu(i)/950))*1.0e+03;
%         C = density .* capacity;
%     elseif (950 < eu(i) && eu(i) <= 1650)
%         density = 4420 - 0.154*(eu(i)-25);
%         density = density * 1e-9; % rescale to [kg/mm^3]
%         capacity = (0.64 + 0.09*((eu(i)-950)/(1650-950)))*1.0e+03;
%         enthalpyJump = 53.0e+03;
%         C =  density .* capacity + ...
%             density .* enthalpyJump * func_pc_der_mart(eu(i), eu_old(i));
%     else
%         density = 3920 - 0.68*(eu(i)-1650);
%         density = density * 1e-9; % rescale to [kg/mm^3]
%         capacity = (0.83)*1.0e+03;
%         enthalpyJump = 286.0e+03;
%         C =  density .* capacity + ...
%             density .* enthalpyJump * func_pc_der_fus(eu(i), eu_old(i));
%     end
% end

% CAPACITY = volumetric heat capacity for SS 316L from ANSYS 2022R2

C = zeros(size(eu));

for i=1:numel(eu)
    if eu(i) <= 1326.9
        density = 7954 - 643*((eu(i)-26.85)/1326.9);
        density = density * 1e-9; % rescale to [kg/mm^3]
        capacity = (498.73 + 172.8*(eu(i)/1326.9));
        C = density .* capacity;
    elseif (1326.9 < eu(i) && eu(i) <= 1426.9)
        density = 7311 - 332*((eu(i)-1326.9)/(1426.9-1326.9));
        density = density * 1e-9; % rescale to [kg/mm^3]
        capacity = 671.53 - 98.33*((eu(i)-1326.9)/(1426.9-1326.9));
        C =  density .* capacity;
    else
        density = 6979 - 848*((eu(i)-1426.9)/(2526.9-1426.9));
        density = density * 1e-9; % rescale to [kg/mm^3]
        capacity = 769.86;
        C =  density .* capacity;
    end
end
end

%Phase Change function and derivative at martensitic-liquid transition
function f_pcDer = func_pc_der_fus(eu, eu_old)

if eu == eu_old
    f_pcDer = 0.0;
else
    f_pcDer = (func_pc_fus(eu) - func_pc_fus(eu_old)) / (eu - eu_old);
end

end

function f_pc = func_pc_fus(eu)

f_pc(eu <= 1650) = 0.0;
f_pc(eu > 1650) = 1.0;

end

%Phase Change function and derivative at solid(alpha)-martensitic(beta) transition
function f_pcDer = func_pc_der_mart(eu, eu_old)

if eu == eu_old
    f_pcDer = 0.0;
else
    f_pcDer = (func_pc_mart(eu) - func_pc_mart(eu_old)) / (eu - eu_old);
end

end

function f_pc = func_pc_mart(eu)

f_pc(eu <= 950) = 0.0;
f_pc(eu > 950) = 1.0;

end

