function out = networkGUI(networkName, link, node, source, sink, path, pathList, nt, t)

% fileName is used when saving a GIF of animation

%% Callback definitions

    function linkSelectionCallback(h,~)
        
        selectedLink = find(hlinks == h);
        flowProfile(selectedLink, link, 'link', t)
        
    end
    function nodeSelectionCallback(h,~)

        selectedNode = find(hnodes == h);
        selectedSource = find(source.nodes == selectedNode);
        selectedSink = find(sink.nodes == selectedNode);
        
        if selectedSource
            flowProfile(selectedSource, source, 'source', t)
        end
        
        if selectedSink
            flowProfile(selectedSink, sink, 'sink', t)
        end
        
    end
    function sliderCallback(h,~)
        
        set(hplay, 'String', 'Start')
        tval = round(get(h, 'Value'));
        updateUi(tval)
        
    end
    function rewindCallback(~,~)
        
        set(hplay, 'String', 'Start')
        simulationInProcess = false;
        tval = 1;
        updateUi(tval)
        
    end
    function popupCallback(h,~)
        
        ind = get(h, 'Value');
        cmaplabel = popupInfo{2,ind};
        cmap = popupInfo{3,ind};
        cind = popupInfo{4,ind};
        caxis(popupInfo{5,ind})
        colormap(cmap)
        hcol = colorbar;
        ylabel(hcol, cmaplabel)
        if strcmp(get(hplay, 'String'), 'Start')
            updateUi(tval)
        end
        
    end
    function displayNodesCallback(h,~)
        
        if get(h, 'Value') == 0
            set(hnodes, 'Visible', 'off')
        else
            set(hnodes, 'Visible', 'on')
        end
        
    end
    function nodeNumbersCallback(h,~)
        
        if get(h, 'Value') == 1
            set(hNodeNumbers, 'Visible', 'on')
        else
            set(hNodeNumbers, 'Visible', 'off')
        end
        
    end
    function linkNumbersCallback(h,~)
        
        if get(h, 'Value') == 1
            set(hLinkNumbers, 'Visible', 'on')
        else
            set(hLinkNumbers, 'Visible', 'off')
        end
        
    end
    function playCallback(h,~)  % *** + GIF generator ***
        
        simulationInProcess = ~simulationInProcess;
        
        if simulationInProcess == true
            set(h, 'String', 'Stop')
            tval = get(hslider, 'Value');
            step = 1;
            
            for tn = tval:step:nt
                if ~simulationInProcess
                    return
                end
                updateUi(tn)
                                %%% start save to GIF %%%
%                                 drawnow
%                                 fileName='DUEsolution';
%                                 frame = getframe(hfig);
%                                 im = frame2im(frame);
%                                 [A,map] = rgb2ind(im,256);
%                                 if tn == 1
%                                     imwrite(A,map,[fileName, '.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
%                                 else
%                                     imwrite(A,map,[fileName, '.gif'],'gif','WriteMode','append','DelayTime',0.2);
%                                 end
                                %%% end save to GIF %%%
                pause(0.1)
            end
            set(h, 'String', 'Start')
            simulationInProcess = false;
            
        else  % if the button string is not 'Start'
            set(h, 'String', 'Start')
        end
        
    end
    function zoomCallback(~,~)
        nscale = get(ax, 'Xlim');
        nscale = nscale(2) - nscale(1);
        zf = (scale / nscale) ^ 0.6;
        set(hlinks, 'LineWidth', zf * linkSize)
        set(hnodes, 'MarkerSize', zf * nodeSize)
        set(hlinks, {'Xdata'}, num2cell(node.X(adjacencyList) + dx(:,[1 1]) / zf, 2), ...
            {'Ydata'}, num2cell(node.Y(adjacencyList) + dy(:,[1 1]) / zf, 2))
        set(hLinkNumbers, {'Position'}, num2cell([xLinkText + dx / zf, ...
            yLinkText + dy / zf, zeros(link.count, 1)], 2))
    end
    function closeReqCallback(~,~)
        
        if simulationInProcess == true
            simulationInProcess = false;
            pause(0.4)
        end
        delete(hfig)
    end
    function pathHighlighterCallback(h,~)
        
        pathNumber = round(str2double(h.String));
        if pathNumber > 0 && pathNumber <= path.count
            pathLinks = pathList(pathNumber, 1:path.linkCount(pathNumber));
            set(hlinks, 'Color', [0.9, 0.9, 1])
            set(hlinks(pathLinks), 'Color', 'b')
            set(hnodes, 'MarkerFaceColor', [0.9, 0.9, 0.9])
            set(hnodes(path.sourceNode(pathNumber)), 'MarkerFaceColor', 'g')
            set(hnodes(path.sinkNode(pathNumber)), 'MarkerFaceColor', 'r')
        else
            h.String = '0';
            set(hlinks, 'Color', [0.9, 0.9, 1])
            set(hnodes, 'MarkerFaceColor', [0.9, 0.9, 0.9])
        end
        
    end
    function pathTravelTimeCallback(~,~)
        
        pathNumber = round(str2double(hpath.String));
        if pathNumber > 0 && pathNumber <= path.count
            fig = figure;
            fig.NumberTitle = 'off';
            fig.Name = ['Path ', num2str(pathNumber), ' Travel Time'];
            ax1 = axes;
            plot(ax1, t / 60, path.travelTime(pathNumber,:) / 60)
            xlabel('Departure time (min)'), ylabel('Travel time (min)')
        end
        
    end

    function updateUi(tn)
        
        set(hslider, 'Value', tn)
        set(hlinks, {'Color'}, num2cell(cmap(cind(:,tn),:), 2))
        set(hnodes, {'MarkerFaceColor'}, num2cell(cmapN(cindN(:,tn),:), 2))
        %         sel = sel1(:,tn)|sel2(:,tn);
        %         set(hnodes(sel), 'MarkerEdge', 'none')
        %         set(hnodes(~sel), 'MarkerEdge', 'k', 'MarkerFace', 'none')
        %         set(hnodes(~sel), 'MarkerFace', 'none')
        textsliderstr = sprintf('%0.1f min', t(tn) / 60);
        set(htextslider, 'String', textsliderstr)
        
    end

%% Create main GUI figure

networkScale = sqrt(node.count);
linkSize = 38 / networkScale;  % link width
nodeSize = 125 / networkScale;  % node diameter

scrsz = get(0,'ScreenSize');
hfig = figure('Name', 'Network Simulation', ...
    'OuterPosition',[5 5 scrsz(3)*0.7 scrsz(4)*0.95], ...
    'NumberTitle', 'off');
hpanelmain = uipanel(hfig, 'Position', [0.2 0 0.8 1], 'Background', 'w');
hpanelside = uipanel(hfig, 'Position', [0 0 0.2 1]);
ax = axes('Parent', hpanelmain, 'Position', [0.05 0.05 0.9 0.9]);
title([networkName, ' Network'], 'FontSize', 22)
hold on

% create virtual nodal coordinates for 2-way flow (left hand drive - for
% right hand drive, reverse the sign for dx and dy.
laneSpacing = (max(node.X) - min(node.X)) / networkScale * 0.05;
phi = atan2(node.Y(link.headNode) - node.Y(link.tailNode), node.X(link.headNode) - node.X(link.tailNode));
dx = sin(phi) * laneSpacing;
dy = -cos(phi) * laneSpacing;

% create coordinates for link numbering
adjacencyList = [link.tailNode, link.headNode];
xLinkText = mean(node.X(adjacencyList), 2);
yLinkText = mean(node.Y(adjacencyList), 2);


%% create colourmap for DTA GUI
% uncongested = [linspace(0.7,0,16)', linspace(1,0.7,16)', linspace(0.8,0,16)'];
% congested = [linspace(1,0.8,48)', linspace(1,0,48)', linspace(0,0,48)'];
% cmaprhor  = [uncongested; congested];
cmaprhor = [linspace(0.7,0,16)', linspace(1,1,16)', linspace(0.8,0,16)'; ...
    linspace(0,1,16)', linspace(1,1,16)', linspace(0,0,16)'; ...
    linspace(1,0.8,32)', linspace(1,0,32)', linspace(0,0,32)'];
uncongested = [linspace(0.7,1,16)', linspace(1,1,16)', linspace(0.8,0,16)'];
congested = [linspace(1,0.8,48)', linspace(1,0,48)', linspace(0,0,48)'];
cmaprho  = [uncongested; congested];
% cmaprho = [linspace(0.2,1,32)',linspace(1,1,32)',linspace(0.2,0,32)'; ...
%     linspace(1,0.8,32)',linspace(1,0,32)',linspace(0,0,32)'];
cmapQ = [linspace(0.9,0,32)',linspace(0.9,0,32)',linspace(0.9,1,32)'; ...
    linspace(0,0,32)',linspace(0,0,32)',linspace(1,0,32)'];
cmapN = colormap([0.8 0.8 0.8; 1 0.2 0.2; 0.2 0.8 0.8; 1 0.2 1]);
cmap = [];

maxrho = max(link.density(:));
maxrhor=max(max(link.density./link.jamDensity(:,ones(1,nt))));
maxQin=max(max(link.Qin./link.capacity(:,ones(1,nt))));
maxQout=max(max(link.Qout./link.capacity(:,ones(1,nt))));
cindrho = min(max(round(size(cmaprho,1) * link.density / maxrho), 1), size(cmaprho,1));
cindrhor = min(max(round(size(cmaprhor,1) * link.density ./ link.jamDensity(:,ones(1,nt))/maxrhor), 1), size(cmaprhor,1)); % >1 & <100
cindQin = max(round(size(cmapQ,1) * link.Qin ./ link.capacity(:,ones(1,nt))/maxQin), 1);  % at least as large as 1
cindQout = max(round(size(cmapQ,1) * link.Qout ./ link.capacity(:,ones(1,nt))/maxQout), 1);
% colour coding for nodes
sel1 = false(node.count,nt); sel2 = sel1;
sel1(source.nodes,:) = source.Qout > 0;
sel2(sink.nodes,:) = sink.Qin > 0;
cindN = ones(size(sel1));
cindN(sel1) = 2; cindN(sel2) = 3; cindN(sel1 & sel2) = 4;
cind = [];

cmaplabelQ = 'Relative Inflow (Q/C %)';
cmaplabelQout = 'Relative Outflow (Q/C %)';
cmaplabelrhor = 'Relative Density (\rho/\rho_jam %)';
cmaplabelrho = 'Density (\rho veh/m)';

ratio = [0 100];
tval = 1;

%% Plot Network
% plot links
hlinks = plot(ax, (node.X(adjacencyList) + dx(:,[1 1]))', ...
    (node.Y(adjacencyList) + dy(:,[1 1]))');
set(hlinks, 'LineWidth', linkSize)

% plot nodes
hnodes = NaN(node.count,1);
for n = 1:node.count
    hnodes(n) = line(node.X(n), node.Y(n));
end

hnodes = sort(hnodes);
set(hnodes, 'Marker', 'o', 'MarkerEdge', 'none', 'HitTest', 'on', 'MarkerFaceColor', [0.2 0.2 0.2])

% display link and node numbers (initially hidden)
hNodeNumbers = text(node.X,node.Y, cellstr(num2str((1:node.count)')), 'horizontal', 'center');
set(hNodeNumbers, 'PickableParts', 'none', 'Visible', 'off')
hLinkNumbers = text(xLinkText+dx, yLinkText+dy, cellstr(num2str((1:link.count)')), 'horizontal', 'center');
set(hLinkNumbers, 'PickableParts', 'none', 'Visible', 'off')

% attach callback functions to nodes and link handles
set(hlinks, 'ButtonDownFcn', @linkSelectionCallback)
set(hnodes, 'ButtonDownFcn', @nodeSelectionCallback)

set(zoom(hfig), 'ActionPostCallback', @zoomCallback, 'Enable', 'on')
scale = get(ax, 'xlim');
scale = scale(2) - scale(1);

%% UI elements
buttonheight = 0.05;
popupInfo(:,1) = {'Absolute Density'; cmaplabelrho; cmapQ; cindrho; [0 maxrho]};
popupInfo(:,2) = {'Relative Density'; cmaplabelrhor; cmapQ; cindrhor; [0 maxrhor*100]};
popupInfo(:,3) = {'Relative flow in'; cmaplabelQ; cmapQ; cindQin; [0 maxQin*100]};
popupInfo(:,4) = {'Relative flow out'; cmaplabelQout; cmapQ; cindQout; [0 maxQout*100]};

% stop animation on close
set(hfig, 'CloseRequestFcn', @closeReqCallback)

% make a slider
hslider = uicontrol(hpanelmain, 'Style', 'slider', ...
    'Position', [100, 5, 500, 20], ...
    'Min', 1, 'Max', nt, ...
    'Value', tval, ...
    'SliderStep', [0.05 0.05], ...
    'Callback', @sliderCallback);

% make text display above slider
htextslider = uicontrol(hpanelmain, 'Style', 'text', ...
    'Background', 'w', ...
    'Position', [5, 7, 90, 20], ...
    'FontSize', 16, ...
    'String', ['t = ', num2str(tval), ' h']);

% make a start/stop button
hplay = uicontrol(hpanelside, 'Style', 'pushbutton', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.85, 0.8, 0.1], ...
    'FontSize', 16, ...
    'String', 'Start', ...
    'Callback', @playCallback);

% make a rewind button
uicontrol(hpanelside, 'Style', 'pushbutton', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.8, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', '|<<', ...
    'Callback', @rewindCallback);

% pop-up choice of density, flow
defpop = 2;
hpop = uicontrol(hpanelside, 'Style', 'popupmenu', ...
    'Units', 'normal', ...
    'Position', [0.05, 0.7, 0.9, buttonheight], ...
    'FontSize', 14, ...
    'String', popupInfo(defpop,:), ...
    'Value', defpop, ...
    'Callback', @popupCallback);

% checkbox for displaying nodes
uicontrol(hpanelside, 'Style', 'checkbox', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.6, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', 'Display nodes', ...
    'Value', 1, ...
    'Callback', @displayNodesCallback);

% checkbox for node numbers
uicontrol(hpanelside, 'Style', 'checkbox', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.5, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', 'Node numbers', ...
    'Callback', @nodeNumbersCallback);

% checkbox for link numbers
uicontrol(hpanelside, 'Style', 'checkbox', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.4, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', 'Link numbers', ...
    'Callback', @linkNumbersCallback);

% Path section
uicontrol(hpanelside, 'Style', 'text', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.3, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', 'Enter path number:');

% path highligh number
hpath = uicontrol(hpanelside, 'Style', 'edit', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.2, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', '0', ...
    'Callback', @pathHighlighterCallback);

% path travel time graph
uicontrol(hpanelside, 'Style', 'pushbutton', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.1, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', 'Path travel time', ...
    'Callback', @pathTravelTimeCallback);


%% initialise UI
popupCallback(hpop,[])
zoomCallback([],[])
axis off equal
simulationInProcess = false;

out = [];

end
