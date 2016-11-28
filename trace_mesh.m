function varargout = trace_mesh(varargin)
% function varargout = trace_mesh(varargin)
%----------------------------------------------------------------
% PURPOSE
%    Plot mesh
%----------------------------------------------------------------
% INPUT
%    mesh structure
%----------------------------------------------------------------
% OUTPUT
%    figure (optional)
%----------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     10-02-2016
%----------------------------------------------------------------
%% Close all initial figures
close all

%% Ouverture des variables globale
global geometry;
global options;

%% Check on varargin
if isempty(varargin)
    
    error('Not enough input arguments');
    
end

if ~mod(nargin,2)==0
    
    error('Input arguments must be odd');
    
end

%% Plotting
h = figure;
box on
hold on
set(gca,'fontsize',14,'FontName','arial')
axis off;

for nb_arg=1:2:(nargin-1)
    
    label = varargin{nb_arg};

    lab_color = varargin{nb_arg+1};
    
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
        % Loop over the elements
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
            
        end
        
        % Figure and related properties
        
        switch options.dimension
            
            case 2
                
                patch(Ex',Ey',lab_color);
                
            case 3
                
                patch(Ex',Ey',Ez',lab_color);
                
        end
        
    end
    
end
legend(varargin{1:2:(nargin-1)},'location','best')
hold off

if options.dimension == 3
    view(3)
end
varargout{1} = h;
