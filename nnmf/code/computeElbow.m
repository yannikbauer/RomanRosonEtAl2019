function [idxOfBestPoint, optVal, fdist] = computeElbow(points,doPlot,xstr,ystr)
% computeElbow(points,xstr,ystr) is a method ofthen used in clutering for
% finding the optimal number of clusters. It computes the optimal datapoint
% in relationship between # of datapoints and datapoints values
%
%   points     = input data
%   doPlot     = displays a figure with the elbow method
%   xstr       = label for the x-axis
%   ystr       = label for the y-axis
%
%   [idxOfBestPoint, optVal, fdist] = computeElbow(...) returns the optimal 
%   value with x = idxOfBestPoint, y = optVal. fdist corresponds to the
%   handles of the figure
%
if nargin < 4
    ystr = '';
end
if nargin < 3
    xstr = '';
end
if nargin < 2
    doPlot = 0;
end
if nargin < 1
    error('computeElbow:Insuficient input!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get coordinates of all the points
nPoints = length(points);
allCoord = [1:nPoints; points']';

%# pull out first point
firstPoint = allCoord(1,:);

% get vector between first and last point - this is the line
lineVec = allCoord(end,:) - firstPoint;

% normalize the line vector
lineVecN = lineVec / sqrt(sum(lineVec.^2));

% find the distance from each point to the line:
% vector between all points and first point
vecFromFirst = bsxfun(@minus, allCoord, firstPoint);

% To calculate the distance to the line, we split vecFromFirst into
% two components, one that is parallel to the line and one that is
% perpendicular. Then, we take the norm of the part that is
% perpendicular to the line and get the distance.
% We find the vector parallel to the line by projecting vecFromFirst
% onto the line. The perpendicular vector is vecFromFirst -
% vecFromFirstParallel. We project vecFromFirst by taking the scalar
% product of the vector with the unit vector that points in the
% direction of the line (this gives us the length of the projection
% of vecFromFirst onto the line). If we multiply the scalar product
% by the unit vector, we have vecFromFirstParallel
scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
vecFromFirstParallel = scalarProduct * lineVecN;
vecToLine = vecFromFirst - vecFromFirstParallel;

% distance to line is the norm of vecToLine
distToLine = sqrt(sum(vecToLine.^2,2));

% plot the distance to the line
% figure('Name','distance from curve to line'), plot(distToLine)

% now all you need is to find the maximum
[~,idxOfBestPoint] = max(distToLine);
optVal = allCoord(idxOfBestPoint,2);

% PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doPlot
    fdist = figure;
    plot(points, 'LineWidth', 1.5);
    hold on
    
    % plot vertical line
    ycor = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);
    xcor = ycor;
    xcor(:) = idxOfBestPoint;
    plot(xcor,ycor,'--', 'color', [0.3 0.3 0.3])
    
    % plot horizontal line
    xcor = linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    ycor = xcor;
    ycor(:) = optVal;
    plot(xcor,ycor,'--', 'color', [0.3 0.3 0.3])
    
    % plot the section with the diagonal; flip the original position
    yVal = min(points):round((max(points) - min(points))/6):max(points);
    
%     set(gca,'XMinorTick','on','YMinorTick','on', 'FontSize', 13)
%     set(gca, 'YTick', round(yVal), 'YTickLabel', round(yVal),'TickDir','out')
%     grid on; grid minor; box 'off'; xlabel(xstr); ylabel(ystr)
%     
%     set(gca, 'ylim', [min(points) max(points)])
%     set(gca, 'xlim', [0 length(points)])
%     
%     title('Distribution')
end
end








