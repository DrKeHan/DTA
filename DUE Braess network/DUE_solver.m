function  [h_old,E,epsilon,iter_needed]=DUE_solver(time_horizon,...
                                                    OD_demand,...
                                                    T_A,...
                                                    ODpath_set,...
                                                    h,...
                                                    Max_iteration,...
                                                    nt,...
                                                    dt,...
                                                    threshold,...
                                                    tolerance,...
                                                    alpha,...
                                                    processedDataFile)
% This function performs the fixed-point algorithm for the simultaneous
% route-and-departure-time choice dynamic user equilibrium model
%
% DOCUMENTATION: 
%        Han,K., Eve,G., Friesz,T.L. (2018) Computing dynamic user equilibria on large-scale
%        networks: From theory to software implementation.
%        arXiv:1810.00777. URL: https://arxiv.org/abs/1810.00777
% 
% *************************
% Developed by Ke Han, 2018
% *************************
                                                
%% Fixed point iteration
epsilon=zeros(1,Max_iteration);
h_old=h;
iter_needed=Max_iteration;

for i=1:Max_iteration
     % network loading 
     fprintf('Starting fixed-point iteration no. %4.0f \n\n', i);
     [Delay]=DYNAMIC_NETWORK_LOADING(h_old, nt, dt, processedDataFile); 
     Delay=Delay/3600; % second to hour;
              
     % fixed point update
     [ E,h_new ] = Update(time_horizon,dt,h_old,Delay,OD_demand,T_A,ODpath_set,alpha,tolerance);

     % check for convergence
     epsilon(i)=sum(sum((h_old-h_new).^2))/(sum(sum(h_old.^2)));
     if epsilon(i)<=threshold
         iter_needed=i;
         epsilon(i+1:end)=[];
         break;
     end
     
     % Adjustment of the fixed point algorithm, see Friesz et al. (2011)
     beta=1/(1+i)^0.9; 
     h_new=beta*h_old+(1-beta)*h_new;  
     h_old=h_new;
end



end
