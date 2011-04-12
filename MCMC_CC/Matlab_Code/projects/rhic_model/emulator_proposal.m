function OUT = emulator_proposal( FLAG, THETA1, THETA2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


global CONFIG;
global RANGES;

F = fieldnames(THETA1);

switch FLAG

case 'rand'
	
	while (true)
		OUT = THETA1;
		for i=1:length(F)
			OUT.(F{i}) = normrnd(OUT.(F{i}),CONFIG.PROPOSAL.MIXING_STDDEV.(F{i}));
		end

		%Checks go here.
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%Recall, some (most) parameters we are optimizing the log...
		NEWPARAMS = OUT;
		for i=1:length(CONFIG.LOGPARAM)
			NEWPARAMS.(CONFIG.LOGPARAM{i}) = exp(NEWPARAMS.(CONFIG.LOGPARAM{i}));
        end

        try
            rangemat = cell2mat(struct2cell(RANGES));
            if (all(arrayfun(@(a,b,c)(a>=b)&&(a<=c),cell2mat(struct2cell(OUT)),rangemat(:,1),rangemat(:,2))==1))
                break;%Seems like everything went OK-- break out of while loop.
            end
            %Otherwise the while loop keeps running...
        catch exception
            OUT
            rangemat(:,1)
            rangemat(:,2)
            exception
            error('Wierd.');
        end
	end

case 'pdf'
	for i=1:length(F)
		OUT(i) = normpdf(THETA1.(F{i}),THETA2.(F{i}),CONFIG.PROPOSAL.MIXING_STDDEV.(F{i}));	
	end
	

otherwise
	error(['Unrecognized flag to proposal',FLAG,'.']);
end


