% This script performs computation of the simultaneous
% route-and-departure-time choice dynamic user equilibrium model
%
% OUTPUTS:
%        path flow vector: h_final
%        path effective delay vector: Eff_delay
%        gap between two consecutive iterates: epsilon
%        number of iterations performed until convergence criterion is met: iter_needed
%        computational time of the DUE solver: elapsedtime
%        
%
% DOCUMENTATION AND CITE AS: 
%        Han,K., Eve,G., Friesz,T.L. (2019) Computing Dynamic User Equilibria on
%        Large-Scale Networks with Software Implementation. Networks and Spatial
%        Economics, https://doi.org/10.1007/s11067-018-9433-y
%        
% 
% *************************
% Developed by Ke Han, 2018
% *************************
 
clc
%% User-defined section
dt=180;           % time step (in second)
threshold=0.001;   % Convergence threshold for the DUE algorithm
tolerence=0;      % Indifference band if considering bounded rationality
Max_iteration=100;  % maximum iteration number for the DUE algorithm
alpha=2;  % time step size of the fixed point algorithm; this measures drivers' sensitivity to the travel costs
processedDataFile = [cd, '/ChicagoSketch250000_pp.mat']; % file path to the pre-processed network file

%% Load relevant data
load(processedDataFile, 'path');
load Network_planning_parameters;  
load OD_info;
Num_OD=size(OD_set,1);  % number of OD pairs
time_horizon=[0, 5]; % time horizon, in hour 
nt=round(range(time_horizon)*3600/dt);  % number of time steps

%% Generating initial path flow vector
%  This section can be commented if supplying user-defined initial path
%  flow vector that satisfies the demand constraints. 

h=zeros(1,nt);  h(round(nt/10):round(nt*2/5))=0.002;
InitialPathFlow=ones(path.count,1)*h;

for i=1:Num_OD
    dummy=InitialPathFlow(ODpath_set{i,1}',:);
    nu=OD_demand(i)/(dt*sum(sum(dummy)));
    dummy=dummy*nu;
    InitialPathFlow(ODpath_set{i,1}',:)=dummy;
end


%% DUE solver
tic;
[h_final,Eff_delay,epsilon,iter_needed]=DUE_solver(time_horizon,... % time horizon
                                           OD_demand,... % OD demand
                                           T_A,... % target arrival times
                                           ODpath_set,... % OD path set
                                           InitialPathFlow,... % initial iterate of the algorithm
                                           Max_iteration,... % maximum iterations of the algorithm
                                           nt,...  % number of time steps
                                           dt,... % time step (in hour)
                                           threshold,... % convergence threshold
                                           tolerence,... % indifference band
                                           alpha,...  % step size
                                           processedDataFile);  % file path to the pre-processed network data
elapsedtime=toc;
OD_gap=zeros(1,Num_OD); % calculate the gap for each OD
for i=1:Num_OD
    P=ODpath_set{i};
    hh=h_final(P,:); EE=Eff_delay(P,:);
    OD_gap(i)=max(EE(hh>=0.5/3600))-min(EE(hh>=0.5/3600)); 
    % the cut-off threshold for the departure rate is 0.5 veh/hr 
end

save('DUE_out.mat', 'dt', 'h_final', 'Eff_delay', 'epsilon', 'iter_needed', 'elapsedtime','OD_gap');
