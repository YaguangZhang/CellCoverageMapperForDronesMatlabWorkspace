classdef AnomalyRegion
    %Represents a geographical region of 1/3 arc-second USGS data
    %   Abstractly represents a matrix of 1x1 degree lat/long "cells".
    %   Various subranges of this elevation data can then be queried at
    %   varying resolutions using the readElevation function.

    properties
        latCellRange % The min/max cell top-left corner latitude range
        longCellRange % The min/max cell top-left corner longitude range
        fltDataFiles % An array of file paths for each cell
    end

    methods
        function [elevData] = readelevation(obj, latRange, longRange, varargin)
        %READELEVATION Reads ElevationData from fltDataFiles
        %   elevData = READELEVATION() reads the entire data range
        %
        %   elevData = READELEVATION(latRange, longRange) reads all values
        %   within the region bounded by latRange and longRange
        %
        %   elevData = READELEVATION(...,'PropertyName',PropertyValue)
        %   specifies various additional properties. The sampleFactor
        %   property specifies how the elevation, latitude, and longitude
        %   matrices/vectors should be sampled when trying to reduce the
        %   data size; i.e., take every nth datapoint (default 1). The
        %   display property is a boolean which specifies whether the debug
        %   print statements should be outputted (default false).
            p = inputParser;
            p.addRequired('obj', @(x) isa(x, 'AnomalyRegion'));
            p.addOptional('latRange', obj.latCellRange - [1, 0], ...
                           @(x) length(x) == 2);
            p.addOptional('longRange', obj.longCellRange + [0, 1], ...
                           @(x) length(x) == 2);
            p.addParameter('sampleFactor', 1, @(x) floor(x) == x && x > 0);
            p.addParameter('display', false, @(x) isa(x, 'logical'));
            p.parse(obj, latRange, longRange, varargin{:});
            inputs = p.Results;

            RASTER_SIZE = 10812;

            [cellLats, cellLongs] = rangeparse(inputs.latRange, ...
                                               inputs.longRange);

            % The appropriate rows and cols to fetch from fltDataFiles
            dataRows = fliplr(inputs.obj.latCellRange(end) - cellLats + 1);
            dataCols = abs(inputs.obj.longCellRange(1)) - abs(cellLongs) + 1;

            sampleSize = ceil(RASTER_SIZE / inputs.sampleFactor);

            elev = zeros(sampleSize * length(cellLats), ...
                         sampleSize * length(cellLongs));
            lats = zeros(sampleSize * length(cellLats), 1);
            longs = zeros(sampleSize * length(cellLongs), 1);

            for i = dataRows
                for j = dataCols
                    if inputs.display
                        fprintf('Reading elevation data %d,%d\n', i, j);
                    end

                    % Read the data file
                    try
                        [e, R] = readgeoraster( ...
                            char(inputs.obj.fltDataFiles(i, j)), ...
                            'CoordinateSystemType', 'geographic');
                    catch
                        error(['Unable to read elevation data. ' ...
                            'Make sure that the specified lat/long ' ...
                            'range lies inside the region.']);
                    end

                    % Sample the elevation data
                    %   e = e(1:inputs.sampleFactor:size(e,1), ...
                    %       1:inputs.sampleFactor:size(e,2));
                    % A workaround for a bug in the terrain elevation
                    % libary. Some downloaded tiles have more data than
                    % what the raster needs. Typically, if the size of the
                    % raster is m x n, the data fetched could be (m+1) x
                    % (n+1). We will only use m x n of the data fetched,
                    % then.
                    e = e(1:inputs.sampleFactor:RASTER_SIZE, ...
                        1:inputs.sampleFactor:RASTER_SIZE);

                    % We only accept GeographicCellsReference for accuracy.
                    if isa(R, 'map.rasterref.GeographicCellsReference')
                        [cellLats, cellLongs] ...
                            = geographicGrid(R, 'gridvectors');
                    else
                        error('Unrecognized raster referencing information.');
                    end

                    % If on the first row of data, copy the cell longitude
                    % values into the overall longitude vector
                    if i == min(dataRows)
                        longs(sampleSize*(j-1) + 1 : sampleSize*j) ...
                            ... = cellLongs(1:inputs.sampleFactor:size(cellLongs,2));
                            = cellLongs(1:inputs.sampleFactor:RASTER_SIZE);
                    end

                    % If on the first column of data, copy the cell
                    % latitude values into the overall latitude vector
                    if j == min(dataCols)
                        lats(sampleSize*(i-1) + 1 : sampleSize*i) ...
                            ... = cellLats(1:inputs.sampleFactor:size(cellLats,2));
                            = cellLats(1:inputs.sampleFactor:RASTER_SIZE);
                    end

                    % Copy the cell elevation matrix into the overall
                    % elevation matrix
                    try
                        elev(sampleSize * (i - 1) + 1 : sampleSize * i, ...
                            sampleSize * (j - 1) + 1 : sampleSize * j) = e;
                    catch err
                        if strcmpi(err.identifier, ...
                                'MATLAB:subsassigndimmismatch')
                            elev(sampleSize * (i - 1) ...
                                + 1 : sampleSize * i, ...
                                sampleSize * (j - 1) ...
                                + 1 : sampleSize * j) ...
                                = e(1:sampleSize, 1:sampleSize);
                        else
                            rethrow(err)
                        end
                    end
                end
            end

            if inputs.display
                fprintf('Trimming data to specified range\n');
            end

            % Trim off latitudes, longitudes, and corresponding elevations
            % that lie outside the specified ranges
            latInds = find(lats > min(inputs.latRange) & ...
                           lats < max(inputs.latRange));
            longInds = find(longs > min(inputs.longRange) & ...
                            longs < max(inputs.longRange));

            elevData = ElevationData;
            elevData.lats = lats(latInds);
            elevData.longs = longs(longInds);
            elevData.elev = elev(latInds, longInds);
        end
    end
end