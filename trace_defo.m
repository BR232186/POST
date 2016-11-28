function varargout = trace_defo(label,lab_color,usol,fact)
% function varargout = trace_defo(label,lab_color,usol,fact)
%--------------------------------------------------------------------------
% PURPOSE
%    Plot deformed mesh
%--------------------------------------------------------------------------
% INPUT
%    label      : label of the geometrical area
%    lab_color  : label of the color
%    usol       : displacement field
%    fact       : amplification factor
%--------------------------------------------------------------------------
% OUTPUT
%    figure (optional)
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     10-02-2016
%--------------------------------------------------------------------------

%% Ouverture des variables globale
global geometry;
global options;
global topology;

%% Check on varargin
if nargin <= 3
    
    error('Wrong number of input arguments');
    
end

%% Plotting
h = figure;
box on
hold on
set(gca,'fontsize',14,'FontName','arial')
axis off;

%% Creation of Ex ... Ez fields
% Number of nodes per element
tmp  = getgeom(label);

if strcmp(tmp.typ_ele,'COMP')
    
    ref_tmp = tmp.ref_comp;
    
    for i=1:length(ref_tmp)
        
        tmp(i) = geometry.objet(ref_tmp(i));
        
    end
    
else
    
    tmp  = getgeom(label);
    
end

for k=1:length(tmp)
    
    % Connectivites
    conn =  tmp(k).elem_conn;
    
    % Initialisation
    Ex = [];
    Ey = [];
    Ez = [];
    
    Ux = [];
    Uy = [];
    Uz = [];
    
    % Loop over the elements
    for j=1:tmp(k).nb_elem
        
        switch options.dimension
            
            case 1
                
                Ux = getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol) * fact;
                
                Ex = [Ex ; geometry.coord(conn(j,:),1)' + Ux'];
                
            case 2

                Ux = getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol) * fact;
                Uy = getdofvaleur_elem(conn(j,:),'RESTR',topology,2,usol) * fact;

                Ex = [Ex ; geometry.coord(conn(j,:),1)' + Ux'];
                Ey = [Ey ; geometry.coord(conn(j,:),2)' + Uy'];
                
            case 3
                
                Ux = getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol) * fact;
                Uy = getdofvaleur_elem(conn(j,:),'RESTR',topology,2,usol) * fact;
                Uz = getdofvaleur_elem(conn(j,:),'RESTR',topology,3,usol) * fact;
                
                Ex = [Ex ; geometry.coord(conn(j,:),1)' + Ux'];
                Ey = [Ey ; geometry.coord(conn(j,:),2)' + Uy'];
                Ez = [Ez ; geometry.coord(conn(j,:),3)' + Uz'];
                
        end
        
    end
    
    % Figure and related properties
    
    switch options.dimension
        
        case 2
            
            patch(Ex',Ey',lab_color);
            
        case 3
            
            patch(Ex',Ey',Ez',lab_color);
            
    end
    
end

legend(label,'location','best')
hold off

if options.dimension == 3
    view(3)
end
varargout{1} = h;
