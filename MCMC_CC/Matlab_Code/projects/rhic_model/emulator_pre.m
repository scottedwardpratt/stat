function [data, rangesstruct] = emulator_pre(ACTUAL)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global CONFIG;

MADAI_ANALYSIS_PATH = '/Users/kevinnovak/Research/RHIC_Research';

setenv('RUN_NAME', 'rhic-2nd-test');

setenv('SCRIPT_HOME', [MADAI_ANALYSIS_PATH '/madai-analysis/src']);

setenv('PATH', [getenv('PATH'), ':/Users/kevinnovak/local:/opt/local/bin:' getenv('SCRIPT_HOME')]);

setenv('DYLD_LIBRARY_PATH', '/Users/kevinnovak/local/lib')

setenv('LD_LIBRARY_PATH', '/Users/kevinnovak/local/lib')

[data, dataerrors] = QueryEmulator([getenv('SCRIPT_HOME') '/../model-data/rhic-2nd-test'], ACTUAL);

filename = [MADAI_ANALYSIS_PATH '/madai-analysis/model-data/' getenv('RUN_NAME') '/ranges.dat'];

[ranges,message] = fopen(filename);

if(ranges == -1)
    filename
    message
    error('Unable to open file.');
end

rangesdata = textscan(ranges, '%*s %s %f %f');

names = rangesdata{1};

for i = 1:size(names,1)
    values{i} = [rangesdata{2}(i) rangesdata{3}(i)];
end
values = values';

rangesstruct = cell2struct(values, names, 1);

% rangesstruct = [rangesdata{2:3}];

end
