% TRANSPARENTIZECURLEGENDS A helper snippet to set the transparency of the
% legend box(es) of the current figure (or subfigures).
%
% Yaguang Zhang, Purdue, 09/13/2017

hLegends = findobj(gcf, 'Type', 'Legend');
for idxLegend = 1:length(hLegends)
    set(hLegends(idxLegend).BoxFace, 'ColorType','truecoloralpha', ...
        'ColorData', uint8(255*[1;1;1;.5]));
end

% EOF