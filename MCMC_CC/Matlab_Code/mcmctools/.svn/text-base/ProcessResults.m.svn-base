function [ Results ] = ProcessResults()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global CONFIG;

%Determine which .mat files to open
numtraces = ceil(CONFIG.OPTIONS.MCMC.MaxT/CONFIG.OPTIONS.MCMC.WriteOut);

tracedata = [];
numdata = 0;

%Copy data in from trace files into a structure array
while(numdata < CONFIG.OPTIONS.STATS.statsamplesize)
    tracefilename = fullfile(CONFIG.OPTIONS.MCMC.OutDir,[CONFIG.META.mcmcrunnickname '_out_' num2str(numtraces) '.mat']);
    current_trace = load(tracefilename,'X');
    numdata = numdata + size(current_trace.X,1);
    tracedata = [tracedata, current_trace]; %concatenating structures forms an array of structures, it doesn't concatenate values
    
    numtraces = numtraces -1;
end
    

F = fieldnames(tracedata(1).X);
results = cell(size(F,1),1);

for param = 1:size(F,1)
    sample = [];
    for i =1:size(tracedata,2)
        currentsample=[tracedata(i).X(:).(F{param})];
        %Form a sample of values where the first value is the last
        %iteration, second value is second to last, etc...
        sample = [sample, fliplr(currentsample)]; 
    end
    sample = sample(1:CONFIG.OPTIONS.STATS.statsamplesize); %...then take only the iterations in our sample.
    %This solves the issue of having a statsample that isn't an integer
    %number of traces
    
    if(any(strcmp(CONFIG.LOGPARAM,F{param}))==1)
        %Exponentiate the parameter values if necessary
        sample = exp(sample);
    end
    
    %Calculate the mean and standard deviation, then store in a cell array
    temp(1)=mean(sample);
    temp(2)=std(sample);
    results{param}=temp;
end

Results = cell2struct(results,F,1);
end

