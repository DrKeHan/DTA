function delay=DYNAMIC_NETWORK_LOADING(pathDepartures, nt, dt, processedDataFile)
% nt is the number of time steps (must be >= number of time steps in the
% departure rate data
% dt is the time step in second

% Notes and assumptions:
% - Link free flow time and kinematic wave travel time are multiples
% of dt (if not in data, then modified before simulating)
% - Sources and sinks represented at OD nodes by means of virtual link
% - All units in SI i.e. metres, seconds and vehicles
%
% DOCUMENTATION AND CITE AS: 
%        Han,K., Eve,G., Friesz,T.L. (2019) Computing Dynamic User Equilibria on
%        Large-Scale Networks with Software Implementation. Networks and Spatial
%        Economics, https://doi.org/10.1007/s11067-018-9433-y
%
% *****************************************
% Developed by Gabriel Eve and Ke Han, 2018
% *****************************************


%% DATA INPUT
load(processedDataFile)
fprintf('Processed data loaded from:      %s\n', processedDataFile)


%% INITIALIZING VARIABLES

% issue warning if any free flow times < time step
if any(link.FFT < dt)
    warning('dt > the minimum link free-flow time\ndt = %1.5f, fft_min = %1.5f', dt, min(link.FFT))
end

% modify array size of pathDepartures to be consistent with nt
nt_departures = size(pathDepartures,2);  % num of time steps in path departure rates
if nt < nt_departures
    nt = nt_departures;            % minimum number of time steps determined by path departure rate data
elseif nt > nt_departures
    pathDepartures(1,nt) = 0;      % pad pathDepartures with zeros to be consistent with nt
end

% round free flow times to multiples of dt and update link lengths accordingly
link.FFT_mod = round(link.FFT / dt) * dt;       % [s] modified free flow time
link.FFT_mod(link.FFT_mod == 0) = dt;           % if free flow time rounded to 0, impose minimum of dt
link.length_mod = link.FFT_mod .* link.length ./ link.FFT;  % [m] modified length of link
N_jam = 4 * link.capacity .* link.FFT_mod;      % (used later in loop)
link.jamDensity = N_jam ./ link.length_mod;     % [veh/m] jam density
tnk = round(link.FFT_mod / dt);                 % num of time-steps for free flow to reach head node
tnw = round(3 * link.FFT_mod / dt);             % num of time-steps for shockwave to propagate to tail node
t = 0:dt:(nt-1)*dt;                             % [s] time vector

Nup = zeros(link.count + source.count + sink.count, nt);    % [veh] cumulative sum of upstream vehicles
Ndn = zeros(link.count + source.count, nt);                 % [veh] cumulative sum of downstream vehicles
Qin = Nup;                                                  % [veh/s] upstream flow
Qout = Ndn;                                                 % [veh/s] downstream flow
Qin_pathLinks = zeros(numTotalPathLinks, nt);               % [veh/s] upstream flow of paths % £££
N_source = zeros(source.count, 1);                          % [veh] number of veh queueing at source at tn


%% Departure Rates

source.departureRates = zeros(source.count,nt);
for n = 1:source.count
    i = source.nodes(n);
    % sum path departure rates for each source
    source.departureRates(n,:) = sum(pathDepartures(path.sourceNode == source.nodes(n), :), 1);
end
% set total departure rate from each source
Qin(source.index,:) = source.departureRates;
Qin_pathLinks(pathSourceLinkIdx,:) = pathDepartures;

clearvars n i

%% COMPUTATION

fprintf('Simulating traffic flows')
dispstat('', 'init')
dispstat(sprintf('\t0%% Complete'))

D_link = zeros(link.count,1);                       % demand = 0 in empty network (initial condition)
D_source = zeros(source.count,1);                   % demand at source (initial)
S_link = link.capacity;                             % supply = capacity (initial condition)

capacity = [link.capacity; inf(source.count,1)];    % flow capacity (including virtual links)
demand = [D_link; D_source];                        % demand (including virtual links)
supply = [S_link; inf(source.count+sink.count,1)];  % supply: how much can enter (including virtual links)

eps = 1e-9;                                         % machine precision tolerance

%tic % timer
for tn = 1:nt  % loop over all time steps
    
    %% Link model
     
    % Link demand
    tkappa = tn - tnk;                  % t - link.length/k index
    ak = find(tkappa >= 1);
    aktk = sub2ind(size(Nup), ak, tkappa(ak));
    
    demandCondition = (Nup(aktk) - Ndn(ak,tn) > eps);  % PDE solution
    
    D_link(:) = 0;
    D_link(ak(demandCondition)) = link.capacity(ak(demandCondition));
    D_link(ak(~demandCondition)) = Qin(aktk(~demandCondition));
    
    % Link supply
    tomega = tn - tnw;                  % t - link.length/w index
    aw = find(tomega >= 1);
    awtw = sub2ind(size(Ndn), aw, tomega(aw));
    
    supplyCondition = (Nup(aw,tn) - Ndn(awtw) - N_jam(aw) < -eps);  % PDE solution
    
    S_link(:) = link.capacity;
    S_link(aw(~supplyCondition)) = Qout(awtw(~supplyCondition));
    

    % Virtual link (source) demand
    isQueue = N_source > eps;               
    D_source(isQueue) = inf;                    % if queue present, demand = infinite
    Q_source = Qin(source.index,tn);            % departure rate
    D_source(~isQueue) = Q_source(~isQueue);    % if no queue, demand = departure rate
    
    
    % Update demand and supply
    demand(link.index) = D_link;
    demand(source.index) = D_source;
    supply(link.index) = S_link;
    
    
    %% Junction model
    
    for i = 1:node.count  % loop through all nodes (junctions)
        nLin = numLinksIn(i);    %  number of incoming links (including the virtual links)
        nLout = numLinksOut(i);  %  number of outgoing links (including the virtual links)
        Lin = linksIn{i};
        
        gamma = cell(nLin,nLout);
        alpha = zeros(nLin,nLout);
        for ik = 1:nLin
            Lik = Lin(ik);
            qin = Qin(Lik,1:tn);
            % logical indexing may be slowing calculation time since not the same size array
            tau = find(qin(Nup(Lik,1:tn) <= Ndn(Lik,tn) ), 1, 'last');
            if isempty(tau)
            % set equal to 1 when there are no values of tau (i.e. before flow has reached end of link)
                tau = 1;
            end
            
            qi = Qin(Lik,tau);
            for jk = 1:nLout
                qijr = Qin_pathLinks(pathLinksIn{i}{ik,jk}, tau);
                % this is because when qi=0 gamma should be 0/is this right?
                if qi<=10^(-10)  % qi is numerically deemed zero, by default all alphas should be zero
                    gamma{ik,jk}=0;
                else
                    gamma{ik,jk} = max(qijr ./ qi, 0);
                    alpha(ik,jk) = sum(gamma{ik,jk});
                end
            end
        end
        
        Lout = linksOut{i};  %  set of outgoing links (including the virtual links) of a given node i
        eta = node.signalPriorities{i};
        effectiveSupply = min([ capacity(Lin), supply(Lout,ones(nLin,1))' ./ alpha ],[], 2);
        qout = min(demand(Lin), eta .* effectiveSupply);
        
        Qin(Lout,tn) = alpha' * qout;           % 	matrix multiplication
        Qout(Lin,tn) = qout;
        
        % split flow back into individual paths
        for ik = 1:nLin
            for jk = 1:nLout
                Qin_pathLinks(pathLinksOut{i}{ik,jk},tn) = gamma{ik,jk} * qout(ik);
            end
        end
        
    end
    
    
    % update number of veh that have entered and exited each link
    Nup(:,tn+1) = Nup(:,tn) + Qin(:,tn) * dt;
    Ndn(:,tn+1) = Ndn(:,tn) + Qout(:,tn) * dt;
    
    % update length of queue at source
    N_source = max(0, N_source + dt*(Qin(source.index,tn) - Qout(source.index,tn)));
    % show percent complete
    if mod(tn, 40) == 0
        dispstat(sprintf('\t%6.0f%% Complete', tn/nt*100))
    end
end

dispstat(sprintf('\t%6.0f%% Complete', tn/nt*100))
%timeElapsed = toc;
%fprintf('\b\tTime elapsed: %6.2fs\n', timeElapsed)

% remove values at Nt+1 for data uniformity
Nup(:,nt+1) = [];
Ndn(:,nt+1) = [];


% clear temp variables
clearvars i n r ak aktk aw awtw nLin nLout Lin Lout eta gamma alpha ik lik qin tau qi


%% Path travel times

fprintf('Computing path travel times\n')

timeFunction = zeros(size(Ndn));
FFT_mod = [link.FFT_mod; zeros(source.count,1)];
% loop through all links, a (virtual links included)
for a = 1:size(Ndn,1)
    % convert free flow times to equivalent number of time-steps
    fft = round(FFT_mod(a)/dt);
    % loop through each time-step, m
    for t_in = 1:nt
        % find link travel time (in # of dt [+1])
        [~,ltt] = min(abs(Ndn(a,t_in:nt) - Nup(a,t_in)));
        if ltt <= fft  % if time of exit is < free flow time
            timeFunction(a,t_in) = t_in + fft;  % travel time = free flow time
        else
            timeFunction(a,t_in) = t_in + ltt-1;  % travel time = travel time
        end
    end
end


path.travelTime = zeros(path.count, nt);
isIncompleteJourney = false;


% loop over all paths
for r = 1:path.count
    % loop over all departure times
    for tn = 1:nt
        % initialise link entry time
        t_in = tn;
        % loop over all links (and virtual link) in path r
        for a = [path.sourceIndex(r), pathList(r,1:path.linkCount(r))]
            % find link exit time / following link entry time
            t_in = timeFunction(a,t_in);
            if t_in > nt
                isIncompleteJourney = true;
                break
            end
        end
        % if vehicle never reaches destination within time horizon
        if isIncompleteJourney
            path.travelTime(r,tn:nt) = (t_in-tn)*dt;
            isIncompleteJourney = false;  % reset condition
            break
        else
            path.travelTime(r,tn) = (t_in - tn) * dt;
        end
    end
end 

clearvars a fft t_in llt isNonArrival r tn


%% OUTPUTS

link.Nup = Nup(link.index,:);
link.Ndn = Ndn(link.index,:);
link.Qin = Qin(link.index,:);
link.Qout = Qout(link.index,:);

source.Nup = Nup(source.index,:);
source.Ndn = Ndn(source.index,:);
source.Qin = Qin(source.index,:);
source.Qout = Qout(source.index,:);

sink.Nup = Nup(sink.index,:);
sink.Qin = Qin(sink.index,:);

% average density on each link
link.density = (link.Nup - link.Ndn) ./ link.length_mod(:, ones(1,nt));


%% Output parameter
delay=path.travelTime;


%% Debugging

% check that total departures is equal~ to total arrivals
totalDepartures = sum(Nup(source.index,end));
totalArrivals = sum(Nup(sink.index,end));
fprintf('Total departures = %0.1f\nTotal arrivals   = %0.1f\n\n', totalDepartures, totalArrivals)





