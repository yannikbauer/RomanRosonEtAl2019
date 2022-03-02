function [figPars, axPars] = aux_setPlotPars()
%   aux_setPlotPars sets parameters for figures and axis
%
%   [figPars, axPars] = aux_setPlotPars() returns the parameters for
%   figures and figure axes

figPars.Color =         'w';
figPars.Units =         'centimeters';
figPars.PaperUnits =    'centimeters';

axPars.Units = 'centimeter';
axPars.xcolor = 'k';
axPars.ycolor = 'k'; 

axPars.TickDir =        'out';
axPars.TickLength =     [0.015 0.025];
axPars.Box =            'off';

axPars.FontSize     = 13;
axPars.LineWidth    = 1.5;
    
end