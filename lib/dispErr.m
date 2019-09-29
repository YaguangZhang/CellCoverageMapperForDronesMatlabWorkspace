function [ ] = dispErr(e)
%DISPERR For displaying catched errors.
%
% Yaguang Zhang, Purdue, 09/18/2019

disp(getReport(e, 'extended', 'hyperlinks', 'on'));

end
% EOF