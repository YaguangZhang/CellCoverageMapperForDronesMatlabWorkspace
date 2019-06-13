function [elevation, resolution] = getElevations(latitude, longitude, varargin)
%GETELEVATIONS queries Google Maps API webservice for ground elevations.
%
%   elevation = GETELEVATIONS(latitude, longitude) returns ground 
%   elevation for latitude and longitude arrays.
%
%   [elevation, resolution] = getElevations(latitude, longitude, ...
%     'key', 'AIzaSyCuN8tjAVEaXorgNjS1tDiVC-oc0QBJoYc' );
%   is an example of a call passing additional attributes to Google Maps 
%   API webservice, and capturing also the resolution of the data. 
%   See https://developers.google.com/maps/documentation/elevation/
%   for details.
%
% EXAMPLE:
%  lat = linspace(36.250278, 36.578581, 40);
%  lon = linspace(-116.825833, -118.291995, 40);
%  dist= linspace(0, 130, 20);
%  elevation = getElevations(lat, lon);
%  plot(dist,elevation)
%  title('Elevation profile from Death Valley to Mt. Whitney in California');
%  xlabel('distance in km');
%  ylabel('elevation in meters');
%
% Author: Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
% License: BSD (Berkeley Software Distribution)
% Documentation: https://developers.google.com/maps/documentation/elevation/

%% process varargin
keyStr = '';
if nargin>2
  p = inputParser;
  p.addParameter('key', '', @ischar)
  p.FunctionName = 'getElevations';
  p.parse(varargin{:})
  results = p.Results;
  if ~isempty(results.key)
    keyStr = sprintf('&key=%s', results.key);
  end
end

%% Check inputs
nPos = numel(latitude);
assert(nPos>0, 'Latitude and longitude inputs can not be empty')
assert(nPos==numel(longitude), 'Latitude and longitude inputs are not of the same length')
assert(min(latitude(:)) >= -90 && max(latitude(:)) <= 90, 'Latitudes has to be between -90 and 90')
assert(min(longitude(:))>=-180 && max(longitude(:))<=180, 'Longitude has to be between -180 and 180')

%% Query Google
elevation  = zeros(size(latitude))*nan;
resolution = zeros(size(latitude))*nan;
batch = [1:50:nPos nPos+1]; % group in batches of 50
for iBatch=2:length(batch)
  idx = batch(iBatch-1):batch(iBatch)-1;
  coord = '';
  for k = 1:length(idx)
    coord = sprintf('%s%9.6f,%9.6f|',coord,latitude(idx(k)),longitude(idx(k)));
  end
  
  %% create query string and run a query
  website = 'https://maps.googleapis.com/maps/api/elevation/xml?locations=';
  url = [website, coord(1:end-1), keyStr];
  str = urlread(url);
  
  %% Parse results
  status = regexp(str, '<status>([^<]*)<\/status>', 'tokens');
  switch status{1}{1}
    case 'OK'
      res = regexp(str, '<elevation>([^<]*)<\/elevation>', 'tokens');
      elevation(idx) = cellfun(@str2double,res);
      if nargout>1
       res = regexp(str, '<resolution>([^<]*)<\/resolution>', 'tokens');
       resolution(idx) = cellfun(@str2double,res);
      end
    case 'INVALID_REQUEST'
      error('Google Maps API request was malformed');
    case 'OVER_QUERY_LIMIT'
      error('Google Maps API requestor has exceeded quota');
    case 'REQUEST_DENIED'
      error('Google Maps API did not complete the request (invalid sensor parameter?)');
    case 'UNKNOWN_ERROR'
      error('Google Maps API: an unknown error.');
  end
end
