function OUT = cosmo_proposal(FLAG,THETA1,THETA2)

global CONFIG;

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
		
		N= CONFIG.STATIC.nz*CONFIG.STATIC.nl; %Total number of data points
		

		V = zeros(1,N);
        %disp('Storing values into V in proposal function.');
		for k = 1:length(CONFIG.CHOL)
    
            if(mod(str2num(CONFIG.CHOL{k}(2)),1)~= 0)
                if(CONFIG.OPTIONS.DEBUG.HandleMessages)
                    HandleMessages(3, 'Error interpreting cholesky parameters.');
                end
                error('Error interpreting cholesky parameters.');
            end
            index = str2num(CONFIG.CHOL{k}(2));
    
            V(index) = NEWPARAMS.(CONFIG.CHOL{k});
        end
        
        
        tempL = [];
        try
            tempL = BuildCholeskyFinal(V);
        catch
            %No need to do anything
            tempL = [];
            if(CONFIG.OPTIONS.DEBUG.HandleMessages)
                HandleMessages(2, 'BuildCholeskyFinal failed. Trying a different proposal set...');
            else
                disp('BuildCholeskyFinal failed.  Trying a different proposal set...');
            end
        end


        if (~isempty(tempL))
            clear tempL;
            break;%Seems like everything went OK-- break out of while loop.
        end
        %Otherwise the while loop keeps running...
	end

case 'pdf'
	for i=1:length(F)
		OUT(i) = normpdf(THETA1.(F{i}),THETA2.(F{i}),CONFIG.PROPOSAL.MIXING_STDDEV.(F{i}));
			
	end
	

otherwise
	error(['Unrecognized flag to proposal',FLAG,'.']);
end
