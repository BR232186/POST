function varargout = reac_curve(tab1,boundary,ev1,idx_var)
% function varargout = reac_curve(tab1,boundary,ev1,idx_var)
%--------------------------------------------------------------------------
% PURPOSE
%    Computation of the reaction curve
%--------------------------------------------------------------------------
% INPUT
%    tab1      : table coming from pasapas
%    boundary  : boundary structure
%    ev1       : loading evolution objects
%    idx_var   : component number to be considered
%--------------------------------------------------------------------------
% OUTPUT
%      h       : figure specifier
%    ev2       : evolution object
%--------------------------------------------------------------------------
% REFERENCES
%     Benjamin RICHARD
%     28-02-2016
%--------------------------------------------------------------------------

%% Ouverture variables
global topology

%% Check on idx_var
if length(idx_var) > 1
   
    error('Only one type of component should be precised')
    
end

%% Initialization
depla = 0;

%% Loop over the time increments
for iincr = 1:length(tab1.resultats)
    
    depla(iincr + 1) = interp1(ev1.valeur_t,ev1.valeur_y,tab1.temps_calc.valeur(iincr + 1));
    
end

%% Re-arrangement of the boundary conditions
charg = gtraite_char(boundary,topology);

%% Loop over the charg structures
idx_dof = [];
for ichar=1:length(charg)
    
    % Looking for the dof number of interest given the local number of the
    % dofs idx_var
    
    switch charg(ichar).type
        
        case 'DEPL'
            
            % Looking for the the dof number
            dimp = charg(ichar).dimp;
            
            % Restriction of the index table
            index_rest = topology.index(dimp(:,1),:);
            
            idx_dof = [idx_dof ...
                index_rest(find(index_rest(:,2) == idx_var),3)];
            
        otherwise
            
            error('Case not impemented');
    
    end
    
end

%% Initialization
react = 0;

%% Loop over the time increments
for iincr = 1:length(tab1.resultats)
    
    react(iincr + 1) = sum(tab1.resultats(iincr).reactions(idx_dof,1));
    
end

%% Plotting
h = figure;
box on
grid on
plot(depla,-react,'b','linewidth',2.0)
xlabel('Displacement (m)','fontsize',18)
ylabel('Load (N)','fontsize',18)
title('Reaction curve','fontsize',18)
box on
grid on
set(gca,'fontsize',18,'fontname','arial')

%% Output storage
l_depla = gmliste('Displacement',depla);
l_react = gmliste('Reaction',react);
ev2     = gmevol('F_delta',l_depla,l_react);

varargout{1} = ev2;
varargout{2} = h;