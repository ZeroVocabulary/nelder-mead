clc;
clear all;

J = @(u) (1 - u(1))^2 + 50*(u(2) - u(1)^2)^2;
J_graph = @(x, y) (1 - x).^2 + 50*(y - x.^2).^2;
% I don't know how to get this number from the function with matlab
% So I made it a part of the arguments
n = 2;
start = [5;5];
printInfo = false;
nm = NelderMead(J, n, start, printInfo);
history = start;
shouldGraph = true;
shouldPause = true;
disp("Press the enter key to see steps");
while(true)
    if(shouldGraph)
        history = showGraph(history, J_graph, nm.x, nm);
    end
    nm.step();
    if(shouldPause), pause; end
    if(nm.shouldTerminate())
        break;
    end
end
disp("Minimum position by manually looping through steps:");
disp(nm.getCurrentMinimum());

nm = NelderMead(J, n, start, printInfo);
disp("Minimum position using minimize function:");
disp(nm.minimize());

disp("Minimum position using fminsearch:");
disp(fminsearch(J,start));

function [history] = showGraph(history, J, p, nm)
    margin = 0.1;
    % Set the range of the graph
    minX = min([p(1,:), history(1,:)]) - margin;
    minY = min([p(2,:), history(2,:)]) - margin;
    maxX = max([p(1,:), history(1,:)]) + margin;
    maxY = max([p(2,:), history(2,:)]) + margin;
    steps = 50;
    x = linspace(minX,maxX,steps);
    y = linspace(minY, maxY,steps);
    [X, Y] = meshgrid(x,y);
    Z = J(X, Y);
    hold off;
    contourf(X,Y,Z,15);
    hold on;
    scatter(p(1,:), p(2,:),5,'r');
    center = calculateCentroid(p);
    history = [history,center];
    plot(history(1,:),history(2,:),'LineWidth',1, 'Color', "yellow");
    minPos = nm.getCurrentMinimum();
    minValue = nm.getValues(minPos);
    title("Min position: " ...
        + round(minPos(1),3) ...
        + ", "...
        + round(minPos(2),3) ...
        + "   Min value: " ...
        + minValue, ...
        "FontSize", ...
        8);
end

function [centroid] = calculateCentroid(p)
    sz = size(p);
    dims = sz(2);
    op = ones([1, dims]) / dims;
    centroid = (op * p')';
end