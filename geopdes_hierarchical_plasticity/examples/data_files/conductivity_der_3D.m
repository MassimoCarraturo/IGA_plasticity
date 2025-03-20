function k_der = conductivity_der_3D(eu)
%CONDUCTIVITY_DER evaluate the temeprature dependent
%conductivity derivative of the material

k_der = cat (1, reshape (zeros(size(eu)), [1, size(eu)]), ...
            reshape (zeros(size(eu)), [1, size(eu)]), ...
            reshape (zeros(size(eu)), [1, size(eu)]));

for i = 1:3
    for j = 1:numel(k_der(i,:))
        if eu(j) > 800
            k_der(i,j) =  34.0;
        else
            k_der(i,j) =  26.7;
        end
    end
end

k_der= k_der*1e-3; %rescale to [W/mm K]

end

