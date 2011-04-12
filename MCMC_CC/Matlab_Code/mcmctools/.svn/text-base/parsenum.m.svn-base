function N = parsenum(N) 
%PARSENUM parse data warehouse mat-file name.
%
%   Example:
%      N = {'out_9_1_10.mat';'out_10_11_1001.mat'};
%      N = repmat(N,100,1);
%

N = strrep(N,'out_',' ');
N = strrep(N,'.mat',';');
N = strrep(N,'_', ',');
N = str2num(cat(2,'[',N{:},']'));
