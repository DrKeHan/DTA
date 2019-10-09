function [ E,h_new,OE ] = Update( time_horizon,dt,h,Delay,demand,T_A,ODpath_set,alpha,tolerance )

% input h is in veh/s
% output h_new is in veh/s
% input dt is in second

h=h*3600;    % convert flow to veh/hour
dt=dt/3600;  % convert dt to hour
time_grid=time_horizon(1):dt:time_horizon(end);time_grid(end)=[];

Num_OD=length(demand);

%% Arrival penalty and effective delay
gamma_early=0.8; gamma_late=1.2;      % coefficients of early and late arrival penalties
AP=zeros(size(Delay));             % initialize arrival penalty
Arrival_Time=ones(size(h,1),1)*time_grid+Delay;
 

for k=1:Num_OD
    for i=1:length(ODpath_set{k,1})
        dummy=Arrival_Time(ODpath_set{k,1}(i),:)-T_A(k);
        dummy(dummy>0)=dummy(dummy>0).^2*gamma_late;
        dummy(dummy<=0)=dummy(dummy<=0).^2*gamma_early;
 
        AP(ODpath_set{k,1}(i),:)=dummy;
    end
end

OE=Delay+AP;   % Original Effective delay (without indifference band)
mu=zeros(1,Num_OD); % minimum effective delay
for i=1:Num_OD
    mu(i)=min(min(OE(ODpath_set{i,1}',:)));
end
 
E=OE; % initialize the revised effective delay E (if indifference band applies)
for i=1:Num_OD
    E(ODpath_set{i,1}',:)=max(E(ODpath_set{i,1}',:),mu(i)+tolerance); % E is the revised effective delay with indifference band
end


%% find dual variable
hh=h-E*alpha;
dualv=zeros(1,Num_OD); 

for i=1:Num_OD
    dualv(i)=fzero(@(zz) userfcn(hh(ODpath_set{i,1}',:),zz,demand(i)),0);
end
 
function [z]=userfcn(x,nu,temp_d)
    temp=x+nu;temp(temp<0)=0;
    z=dt*sum(sum(temp))-temp_d;
end
 
%%
h_new=zeros(size(h));
for i=1:Num_OD
    h_new(ODpath_set{i,1}',:)=h(ODpath_set{i,1}',:)-E(ODpath_set{i,1}',:)*alpha+dualv(i);
end

h_new(h_new<0)=0;
h_new=h_new/3600; % convert flow from veh/hour to veh/s
 
end
 
 

