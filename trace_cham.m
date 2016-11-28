function varargout = trace_cham(idx,mo,Cham_in,field_name,idx_var,varargin)
% function varargout = trace_cham(idx,mo,Cham_in,field_name,idx_var,varargin)
%--------------------------------------------------------------------------
% PURPOSE
%    Plot isovalue on a mesh from a element wise defined field
%--------------------------------------------------------------------------
% INPUT
%    idx        : keyword
%                   - INITIAL  : isovalue in the initial mesh
%    mo         : model structure
%    Cham_in    : Cham structure
%    field_name : name of the field to be plotted
%           . eps0
%           . var0
%           . sig0
%           . epsf
%           . varf
%           . sigf
%    idx_var    : local index of the component to be plotted (depending on
%                 the dimension and on the number of internal variables if
%                 applicable)
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

%% Plotting
h = figure;
box on
hold on
set(gca,'fontsize',14,'FontName','arial')
axis off;

%% Extrapolation of the field
Field_out = extrap_elem_node(mo,Cham_in,field_name,idx_var);

%% Creation of Ex ... Ez fields

for imod=1:length(mo)
    
    % Number of nodes per element
    tmp  = getgeom(mo(imod).lieu);
    
    
    if strcmp(tmp.typ_ele,'COMP')
        
        ref_tmp = tmp.ref_comp;
        
        for i=1:length(ref_tmp)
            
            tmp(i) = geometry.objet(ref_tmp(i));
            
        end
        
    else
        
        tmp  = getgeom(mo(imod).lieu);
        
    end
    
    for k=1:length(tmp)
        
        % Connectivites
        conn =  tmp(k).elem_conn;
        
        % Initialisation
        Ex = [];
        Ey = [];
        Ez = [];
        
        Ux = [];
        
        % Loop over the elements
        switch idx
            
            case 'INITIAL'
                
                for j=1:tmp(k).nb_elem
                    
                    switch options.dimension
                        
                        case 1
                            
                            Ex = [Ex ; geometry.coord(conn(j,:),1)'];
                            
                        case 2
                            
                            Ex = [Ex ; geometry.coord(conn(j,:),1)'];
                            Ey = [Ey ; geometry.coord(conn(j,:),2)'];
                            
                        case 3
                            
                            Ex = [Ex ; geometry.coord(conn(j,:),1)'];
                            Ey = [Ey ; geometry.coord(conn(j,:),2)'];
                            Ez = [Ez ; geometry.coord(conn(j,:),3)'];
                            
                    end
                    
                    Ux = [Ux ; Field_out(conn(j,:),1)'];
                    
                end
                
            otherwise
                
                error('Case not implemented')
                
        end
        
        % Figure and related properties
        switch options.dimension
            
            case 2
                
                patch(Ex',Ey',Ux');
                
            case 3
                
                patch(Ex',Ey',Ez',Ux');
                
        end
        
    end
    
end

colorbar
hold off

if options.dimension == 3
    view(3)
end
varargout{1} = h;
