function f = moving_heat_source_3D(P,laser_radius, laser_depth)

%Gaussian volumentric heat source with laser power P
f = @(x,y,z,path_x,path_y,path_z) 6.0*sqrt(3.0)*P/(pi * sqrt(pi) * laser_radius^2 * laser_depth) * ...
    exp(- 3*(x-path_x).^2/laser_radius^2 - 3*(y-path_y).^2/laser_radius^2 - 3*(z-path_z).^2/laser_depth^2);

% f = @(x,y,z,path_x,path_y,path_z,t) 3.0*5830.0/(pi * 0.015 * 0.010) * ...
%     exp(- 3*(x-path_x).^2/0.015^2 - 3*(y-path_y).^2/0.010^2);

end


