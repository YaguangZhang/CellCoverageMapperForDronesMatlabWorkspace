%TESTITM Test built-in ITM model implementations.
%
% Yaguang Zhang, Purdue, 08/08/2022

pmItm = propagationModel('longley-rice', ...
    'TimeVariabilityTolerance', 0.95, ...
    'SituationVariabilityTolerance', 0.95);

% There is no clear way of feeding in terrain profiles for point-to-point
% links.
disp(pmItm);

% EOF