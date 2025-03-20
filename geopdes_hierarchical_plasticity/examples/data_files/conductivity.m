function k = conductivity(eu)
%CONDUCTIVITY evaluate the temperature dependent
%conductivity of the material

k = zeros(size(eu));

for i=1:numel(eu)
    if eu(i) <= 1370.9
        k(i) =  12.97 + 19.44 * (eu(i)/1370.9);
    elseif (1370.9 < eu(i) && eu(i) <= 1426)
        k(i) =  32.41 - 5.17 * ((eu(i)-1370.9)/(1426-1370.9));
    else
        k(i) =  27.24;
    end
end

% for i=1:size(eu,1)
%     for j=1:size(eu,2)
%         if eu(i,j) >= 800
%             k(i,j) =  34.0 * (eu(i,j)/800);
%         else
%             k(i,j) =  34; % 26.7 * ((800-eu(i,j))/800);
%         end
%     end
% end

k= k*1e-3; %rescale to [W/mm K]

end

