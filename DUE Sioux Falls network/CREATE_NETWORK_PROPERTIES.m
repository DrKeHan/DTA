% CREATE NETWORK PROPERTIES [1]
% - Script for bundling network data into a MAT file in specified
% format for use with the dynamic network loading (DNL) program

clear
clear

%% DATA INPUT

% SOURCE FLOW PRIORITY
sourceSignalPriority = 0.1;                % TODO: arbitrary value

% PROVIDE FILEPATHS FOR NETWORK AND PATH DATA HERE
networkData = [cd, '/SiouxFalls_dat.mat'];
pathData = [cd, '/SiouxFalls6180_paths.mat'];


% load variables from specified MAT-files
load(networkData, 'networkName', 'linkData', '*')
load(pathData, 'pathList')
fprintf('Network data loaded from: %s\n', networkData)
fprintf('Path data loaded from:    %s\n', pathData)

%% Adjacency list/matrix

adjacencyList = linkData(:,1:2);
link.tailNode = adjacencyList(:,1);
link.headNode = adjacencyList(:,2);
link.count = size(linkData, 1);
node.count = length(unique(adjacencyList(:)));
if node.count ~= max(adjacencyList(:));
    warning('Nodes either do not start at ''1'' or are not sequential')
end

% obtain matrix from adjacency list
adjacencyMatrix = sparse(link.tailNode, link.headNode, true, node.count, node.count);


%% Link properties

link.capacity = linkData(:,3);                  % [veh/s] link capacity
link.length = linkData(:,4);                  % [m] original length of link
link.FFT = linkData(:,5);                % [s] original free flow time


%% Nodal coordinates (* if available)

if exist('nodeCoordinates', 'var')
    node.X = nodeCoordinates(:,1);
    node.Y = nodeCoordinates(:,2);
end

%% Signal priorities (+ fixed-alphas) $$$

% create sets of ingoing and outgoing links for each node and signal
% priorities based on ratio of link capacities

[~,c] = find(pathList');  % list path index for each path item % £££
c_ = find([diff(c); 1]);  % find index at end of each path % £££
path.linkCount = [c_(1); diff(c_)];  % find number of links in each path % £££

path.sourceNode = link.tailNode(pathList(:,1)); % £££
path.sinkNode = link.headNode(pathList(sub2ind(size(pathList),1:size(pathList,1),path.linkCount'))); % £££
[source.nodes, ~, pathSources] = unique(path.sourceNode); % £££
[sink.nodes, ~, pathSinks] = unique(path.sinkNode); % £££
% source.nodes = unique(link.tailNode);  % $$$
% sink.nodes = unique(link.headNode); % $$$

source.count = numel(source.nodes);
sink.count = numel(sink.nodes);
path.count = size(pathList,1);

link.index = 1:link.count;
source.index = link.count + (1:source.count);
sink.index = source.index(end) + (1:sink.count);
path.sourceIndex = source.index(pathSources);
path.sinkIndex = sink.index(pathSinks);

node.signalPriorities = cell(node.count,1);
% alphas = node.signalPriorities;                 % $$$  % pre set turning ratios
linksIn = cell(node.count,1);             % list of links into each node (including sources)
linksOut = cell(node.count,1);            % list of links out of each node (including sinks)
numLinksIn = zeros(node.count,1);         % number of flows into node
numLinksOut = zeros(node.count,1);        % number of flows out of node
linkInIdx = zeros(source.index(end),1);
linkOutIdx = zeros(sink.index(end),1);

B = [link.headNode; source.nodes];
A = [link.tailNode; zeros(size(source.nodes)); sink.nodes];
for i = 1:node.count
    Lin = find(B == i);
    Lout = find(A == i);
    linksIn{i} = Lin;
    linksOut{i} = Lout;
    
    nLin = numel(Lin);
    nLout = numel(Lout);
    numLinksIn(i) = nLin;
    numLinksOut(i) = nLout;
    
    linkInIdx(Lin) = 1:nLin;
    linkOutIdx(Lout) = 1:nLout;
    
    % % fixed alphas
    % % simple alpha generator (all equal for every junction) % $$$
    %     alphas{i} = repmat(1/nLout, nLin, nLout);
    
    if sum(source.nodes == i) == 1            % if node i is a source
        Cin = link.capacity(Lin(1:nLin-1));
        node.signalPriorities{i} = [Cin / sum(Cin) * (1-sourceSignalPriority); sourceSignalPriority];
        
    else
        Cin = link.capacity(Lin);
        node.signalPriorities{i} = Cin / sum(Cin);
    end
end

% clear temp variables
clearvars A B Lin Lout nLin nLout Cin

node2source = zeros(node.count,1);
node2sink = zeros(node.count,1);
node2source(source.nodes) = 1:source.count;  % sparse?
node2sink(sink.nodes) = 1:sink.count;

%% Paths  % £££

% create individual index for path flow on every link
% pre-allocate
numPathLinks = numel(find(pathList));
pathSourceLinkIdx = numPathLinks + (1:path.count);
pathSinkLinkIdx = pathSourceLinkIdx(end) + (1:path.count);
ipl = 0;  % initialise pathlink index counter

sourceInIdx = linkInIdx(path.sourceIndex); % link-in index of each path source
sinkOutIdx = linkOutIdx(path.sinkIndex);  % link-out index of each path sink

% process paths
pathLinksIn = cell(node.count,1);
for i = 1:node.count
    pathLinksIn{i} = cell(numLinksIn(i), numLinksOut(i));
end
pathLinksOut = pathLinksIn;

for r = 1:path.count
    pl = path.linkCount(r);
    pth = pathList(r,1:pl);
    kin = sourceInIdx(r);
    kout = linkOutIdx(pth(1));
    pathLinksIn{path.sourceNode(r)}{kin,kout}(end+1) = pathSourceLinkIdx(r);
    ipl = ipl + 1;
    pathLinksOut{path.sourceNode(r)}{kin,kout}(end+1) = ipl;
    for i = 1:pl-1
        Lin = pth(i);
        Lout = pth(i+1);
        kin = linkInIdx(Lin);
        kout = linkOutIdx(Lout);
        pathLinksIn{link.headNode(Lin)}{kin,kout}(end+1) = ipl;
        ipl = ipl + 1;
        pathLinksOut{link.headNode(Lin)}{kin,kout}(end+1) = ipl;
    end
    kin = linkInIdx(pth(pl));
    kout = sinkOutIdx(r);
    pathLinksIn{path.sinkNode(r)}{kin,kout}(end+1) = ipl;
    pathLinksOut{path.sinkNode(r)}{kin,kout}(end+1) = pathSinkLinkIdx(r);
end

numTotalPathLinks = pathSinkLinkIdx(end);

% clear temp variables
clearvars c c_ i ipl kin kout Lin Lout pth pl r sourceInIdx sinkInIdx



%% Save to file

% create fileName string without space chars
fileName = networkName;
fileName(fileName == ' ') = '';  % remove space chars

% save pre processed data with name format 'network_numOfPaths_pp.mat'
dir = [cd, '/', fileName, num2str(path.count), '_pp.mat'];
save(dir)
fprintf(['Variables saved to:       ', dir, '\n\n'])


