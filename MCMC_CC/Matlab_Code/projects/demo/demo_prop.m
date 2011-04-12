function [ OUT ] = demo_prop( FLAG, THETA1, THETA2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global CONFIG;

F=fieldnames(THETA1);

switch FLAG
    case 'rand'
        OUT = THETA1;
        for i=1:length(F)
            STDDEV=CONFIG.PROPOSAL.MIXING_STDDEV(F{i});
            OUT.(F{i})=normrnd(OUT.(F{i}),STDDEV);
        end
    case 'pdf'
    otherwise
        error(['Unrecognized flag to proposal',FLAG,'.']);
end
end

