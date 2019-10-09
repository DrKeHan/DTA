% openNetworkGUI

if isfield(node, 'X')
    
    dispstat('', 'init')
    dispstat('Opening GUI...')
    s = networkGUI(networkName, link, node, source, sink, path, pathList, nt, t);
    dispstat(sprintf('GUI Loaded\n'))
else
    fprintf('Missing node coordinates\n')
end

