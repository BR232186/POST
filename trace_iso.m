function varargout = trace_iso(idx,label,usol,type,varargin)
% function varargout = trace_iso(idx,label,usol,type,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Plot isovalue of a nodal field on a mesh
%--------------------------------------------------------------------------
% INPUT
%    idx        : keyword
%                   - DEFORME  : isovalue on the deformed mesh
%                   - INITIAL  : isovalue in the initial mesh
%    label      : label of the geometrical area
%    usol       : displacement field
%    type       : local type of degrees of freedom to be plotted
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
if nargin < 3
    
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
    switch idx
        
        case 'INITIAL'
            
            for j=1:tmp(k).nb_elem
                
                switch options.dimension
                    
                    case 1
                        
                        Ux = [Ux ; getdofvaleur_elem(conn(j,:),'RESTR',1,usol)'];
                        
                        Ex = [Ex ; geometry.coord(conn(j,:),1)'];
                        
                    case 2
                        
                        Ux = [Ux ; getdofvaleur_elem(conn(j,:),'RESTR',1,usol)'];
                        Uy = [Uy ; getdofvaleur_elem(conn(j,:),'RESTR',2,usol)'];
                        
                        Ex = [Ex ; geometry.coord(conn(j,:),1)'];
                        Ey = [Ey ; geometry.coord(conn(j,:),2)'];
                        
                    case 3
                        
                        Ux = [Ux ; getdofvaleur_elem(conn(j,:),'RESTR',1,usol)'];
                        Uy = [Uy ; getdofvaleur_elem(conn(j,:),'RESTR',2,usol)'];
                        Uz = [Uz ; getdofvaleur_elem(conn(j,:),'RESTR',3,usol)'];
                        
                        Ex = [Ex ; geometry.coord(conn(j,:),1)'];
                        Ey = [Ey ; geometry.coord(conn(j,:),2)'];
                        Ez = [Ez ; geometry.coord(conn(j,:),3)'];
                        
                end
                
            end
            
        case 'DEFORME'
            
            fact = varargin{1};
            
            for j=1:tmp(k).nb_elem
                
                switch options.dimension
                    
                    case 1
                        
                        Ux = [Ux ; getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol)' * fact];
                        
                        Ex = [Ex ; geometry.coord(conn(j,:),1)' + ...
                            getdofvaleur_elem(conn(j,:),'RESTR',1,usol)' * fact];
                        
                    case 2
                        
                        Ux = [Ux ; getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol)' * fact];
                        Uy = [Uy ; getdofvaleur_elem(conn(j,:),'RESTR',topology,2,usol)' * fact];
                        
                        Ex = [Ex ; geometry.coord(conn(j,:),1)' + ...
                            getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol)' * fact];
                        Ey = [Ey ; geometry.coord(conn(j,:),2)' + ...
                            getdofvaleur_elem(conn(j,:),'RESTR',topology,2,usol)' * fact];
                        
                    case 3
                        
                        Ux = [Ux ; getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol)' * fact];
                        Uy = [Uy ; getdofvaleur_elem(conn(j,:),'RESTR',topology,2,usol)' * fact];
                        Uz = [Uz ; getdofvaleur_elem(conn(j,:),'RESTR',topology,3,usol)' * fact];
                        
                        Ex = [Ex ; geometry.coord(conn(j,:),1)' + ...
                            getdofvaleur_elem(conn(j,:),'RESTR',topology,1,usol)' * fact];
                        Ey = [Ey ; geometry.coord(conn(j,:),2)' + ...
                            getdofvaleur_elem(conn(j,:),'RESTR',topology,2,usol)' * fact];
                        Ez = [Ez ; geometry.coord(conn(j,:),3)' + ...
                            getdofvaleur_elem(conn(j,:),'RESTR',topology,3,usol)' * fact];
                        
                end
                
            end
            
    end
    
    % Figure and related properties
    switch options.dimension
        
        case 2
            
            switch type
                
                case 1
                    patch(Ex',Ey',Ux' / fact);
                case 2
                    patch(Ex',Ey',Uy' / fact);
            end
            
            
        case 3
            switch type
                
                case 1
                    patch(Ex',Ey',Ez',Ux' / fact);
                case 2
                    patch(Ex',Ey',Ez',Uy' / fact);
                case 3
                    patch(Ex',Ey',Ez',Uz' / fact);
            end
            
    end
    
end

colorbar
hold off

if options.dimension == 3
    view(3)
end
varargout{1} = h;
