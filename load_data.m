function [ym, tm, ya, ta] = load_data(file)
% Load measurement data from 'Sensor Fusion' smartphone app
%
% USAGE
%   ym = load_data(file)
%   [ym, tm, ya, ta] = load_data(file)
%
% DESCRIPTION
%   Loads magnetometer and accelerometer measurement data collected using
%   the 'Sensor Fusion' smartphone app from the specified file. The
%   magnetometer data is returned in ym and the magnetometer timestamps in
%   tm, the accelerometer data in ya and their timestamps in ta.
%
% PARAMETERS
%   file    Filename.
%
% RETURNS
%   ym      3xNm matrix of magnetometer measurement data. The rows are the
%           three magnetometer axes (x, y, z), the columns the samples in 
%           time.
%   tm      1xNm vector of timestamps for the magnetometer measurements.
%   ya      3xNa matrix of accelerometer measurement data. The rows are the
%           three accelerometer axes (x, y, z), the columns the samples in
%           time.
%   ta      1xNa vector of timestamps for the accelerometer measurements.
%
% AUTHOR
%   2018-08-29 -- Roland Hostettler <roland.hostettler@aalto.fi>

    % Open data file and read all data
    narginchk(1, 1);
    fp = fopen(file, 'r');
    if fp ~= -1
        %data = textscan(fp, '%u64%s%f%f%f', 'delimiter', '\t');
        data = textscan(fp, '%f%s%f%f%f%f%f%f', 'delimiter', '\t');
        fclose(fp);
    else
        error('Can not open file %s for reading.', file);
    end

    % Convert timestamps to seconds, offset at zero
    t = data{1}.'*1e-3;
    t = t-min(t);
    y = cat(2, data{3:5}).';

    % Get accelerometer data
    ia = strcmp(data{2}, 'ACC');
    ta = t(ia);
    ya = y(:, ia);

    % Get magnetometer data
    im = strcmp(data{2}, 'RAWMAG');
    tm = t(im);
    ym = y(:, im);
end