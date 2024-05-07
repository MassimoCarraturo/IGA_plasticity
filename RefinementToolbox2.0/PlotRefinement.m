function fig = PlotRefinement(basis,hspace,X,Y,Z,hpoints)

assert(basis.dim == 2,'Only dimesion of 2 is implemented')

func = ones(size(X));

% Plot the surface

fig = figure;
if exist('hpoints','var')
    surf(X,Y,Z,func,'edgecolor','none')
    colormap summer
    alpha 0.5
    axis equal
    hold on
    points = 1;
else
    % Plot the surface
    surf(X,Y,Z,func,'edgecolor','none')
    colormap summer
    axis equal
    hold on
    points = 0;
end

% Plot geometry boundary
xedge = X(1,:);
yedge = Y(1,:);
zedge = Z(1,:);
plot3(xedge,yedge,zedge,'k')

xedge = X(:,end)';
yedge = Y(:,end)';
zedge = Z(:,end)';
plot3(xedge,yedge,zedge,'k')

xedge = X(:,1)';
yedge = Y(:,1)';
zedge = Z(:,1)';
plot3(xedge,yedge,zedge,'k')

xedge = X(end,:);
yedge = Y(end,:);
zedge = Z(end,:);
plot3(xedge,yedge,zedge,'k')


% Plot all initial elements
resolution = size(X)-1;

nu = resolution(1)/basis.nelem(1);
nv = resolution(2)/basis.nelem(2);

for e2 = 1:basis.nelem(2)
    for e1 = 1:basis.nelem(1)
               
        % Corner points
        u1 = (e1-1)*nu+1;
        u2 = e1*nu+1;
        v1 = (e2-1)*nv+1;
        v2 = e2*nv+1;
        
        xelem = X(u1,v1:v2);
        yelem = Y(u1,v1:v2);
        zelem = Z(u1,v1:v2);
        plot3(xelem,yelem,zelem,'k')
        
        xelem = X(u1:u2,v2)';
        yelem = Y(u1:u2,v2)';
        zelem = Z(u1:u2,v2)';
        plot3(xelem,yelem,zelem,'k')
        
        xelem = X(u2:u1,v2)';
        yelem = Y(u2:u1,v2)';
        zelem = Z(u2:u1,v2)';
        plot3(xelem,yelem,zelem,'k')
        
        xelem = X(u2,v2:v1);
        yelem = Y(u2,v2:v1);
        zelem = Z(u2,v2:v1);
        plot3(xelem,yelem,zelem,'k')
    end
end


% Define refined element structure
hmesh = MakeHmesh(basis,hspace);

% Plot the finer level elements
for L = 1:hspace.level-1
    nelems = basis.nelem*2^(L-1);
    n = 0;
    nu = resolution(1)/nelems(1);
    nv = resolution(2)/nelems(2);
    
    for e2 = 1:nelems(2)
        for e1 = 1:nelems(1)
            n = n+1;
            if hmesh{L}(n) == 1
                
                % Corner points
                u1 = (e1-1)*nu+1;
                u2 = e1*nu+1;
                v1 = (e2-1)*nv+1;
                v2 = e2*nv+1;
                % Middle points
                um = round((u1+u2)/2);
                vm = round((v1+v2)/2);
                
                xsplit = X(um,v1:v2);
                ysplit = Y(um,v1:v2);
                zsplit = Z(um,v1:v2);
                plot3(xsplit,ysplit,zsplit,'k')
                
                xsplit = X(u1:u2,vm);
                ysplit = Y(u1:u2,vm);
                zsplit = Z(u1:u2,vm);
                plot3(xsplit,ysplit,zsplit,'k')
            end

        end
    end
end



% plot control points if data is given

if points  
    
    assert(numel(hspace.ref) == numel(hpoints),'Error in amount of refinements in hspace.ref and hpoints')
    
    nref = numel(hspace.ref);
    cpts = hpoints{nref};
    dim  = size(hpoints{nref},2);

    assert(dim <= 3, 'Dimension higher than 3 cannot be plotted')

    if dim == 1
        scatter(cpts,zeros(1,length(cpts)),'b')
    elseif dim == 2
        scatter(cpts(:,1),cpts(:,2),'b')
    elseif dim == 3
        scatter3(cpts(:,1),cpts(:,2),cpts(:,3),'b')
    end

    if nref > 1
        n = 0;
        iref = [];
        for L = 1:numel(hspace.ref{nref-1}.active)

            % newL is 1 for all new functions in that level
            funcsnref    = hspace.ref{nref}.active{L} + hspace.ref{nref}.deactivated{L};
            funcsnrefmin = hspace.ref{nref-1}.active{L} + hspace.ref{nref-1}.deactivated{L};
            if any(funcsnref - funcsnrefmin)
                funcsminus = funcsnref-2*funcsnrefmin;
                funcsonly  = funcsminus(find(funcsminus));

                iref = [iref find(funcsonly+1)+n];

            end
            n = n + length(find(hspace.ref{nref}.active{L}));
        end
        if numel(hspace.ref{nref}.active) > numel(hspace.ref{nref-1}.active)
            iref = [iref (1:length(find(hspace.ref{nref}.active{L+1})))+n]; 
        end

        if dim == 1
            scatter(cpts(iref),zeros(1,length(cpts(iref))),'r')
        elseif dim == 2
            scatter(cpts(iref,1),cpts(iref,2),'r')
        elseif dim == 3
            scatter3(cpts(iref,1),cpts(iref,2),cpts(iref,3),'r')
        end
    end
    
end

end
                
%                 uind = (e1-1)*nu+1:e1*nu;
%                 vind = (e2-1)*nv+1:e2*nv;

%                 xsplit = [X(round((e2-1)*nv+.5*nv),uind)';X(vind,round((e1-1)*nu+.5*nu))];
%                 ysplit = [Y(round((e2-1)*nv+.5*nv),uind)';Y(vind,round((e1-1)*nu+.5*nu))];
%                 zsplit = [Z(round((e2-1)*nv+.5*nv),uind)';Z(vind,round((e1-1)*nu+.5*nu))];
%                 plot3(xsplit,ysplit,zsplit,'k')