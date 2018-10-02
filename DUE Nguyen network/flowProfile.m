function flowProfile(selectedObject, object, type, t)

% plot Nup, Ndn, Qin, Qout

% units
sec = 1;
secSymbol = 'Time (s)';
hour = 3600;            % 1 hour = 3600 seconds
hourSymbol = 'Time (h)';

tNorm = t / hour;
timeUnitText = hourSymbol;

fig = figure;
fig.NumberTitle = 'off';
ax1 = subplot(2,1,1);
hold on, box off,
ax1.YGrid = 'on';
xlabel(timeUnitText), ylabel('No. of veh')
ax2 = subplot(2,1,2);
hold on, box off
xlabel(timeUnitText), ylabel('Flow (veh/s)')

switch type
    case 'link'
        fig.Name = ['Link ', num2str(selectedObject)];
        plot(ax1, tNorm, object.Nup(selectedObject, :), 'b')
        plot(ax1, tNorm, object.Ndn(selectedObject, :), 'r')
        legend(ax1, 'N_{up}', 'N_{down}', 'Location', 'Best')
        
        plot(ax2, tNorm, object.Qin(selectedObject, :), 'b')
        plot(ax2, tNorm, object.Qout(selectedObject, :), 'r')
        legend(ax2, 'Q_{in}', 'Q_{out}', 'Location', 'Best')
        
    case 'source'
        fig.Name = ['Source ', num2str(selectedObject), ...
            ' (Node ', num2str(object.nodes(selectedObject)), ')'];
        plot(ax1, tNorm, object.Nup(selectedObject, :), 'b')
        plot(ax1, tNorm, object.Ndn(selectedObject, :), 'r')
        legend(ax1, 'N_{up}', 'N_{down}', 'Location', 'Best')
        
        plot(ax2, tNorm, object.Qin(selectedObject, :), 'b')
        plot(ax2, tNorm, object.Qout(selectedObject, :), 'r')
        legend(ax2, 'Q_{in}', 'Q_{out}', 'Location', 'Best')
        
    case 'sink'
        fig.Name = ['Sink ', num2str(selectedObject), ...
            ' (Node ', num2str(object.nodes(selectedObject)), ')'];
        plot(ax1, tNorm, object.Nup(selectedObject, :), 'b')
        
        plot(ax2, tNorm, object.Qin(selectedObject, :), 'b')
        
    otherwise
        error('Invalid type. Accepted types: link, source, sink')
end
