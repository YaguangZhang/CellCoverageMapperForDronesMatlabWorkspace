function hOut = autoText(textString,varargin)
% Function like text() but with auto placement like legend:
%    https://www.mathworks.com/matlabcentral/answers/338565-function-like-text-but-with-auto-placement-like-legend
l = legend(textString,varargin{:});
t = annotation('textbox');
t.String = textString;
t.Position = l.Position;
delete(l);
t.LineStyle = 'None';
if nargout
    hOut = t;
end
% EOF