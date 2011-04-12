function [X,RATIO,AUG] = mcmc(Xinit, varargin)

%MCMC Markov Chain Monte Carlo parameter estimation.
%
%Theta = MCMC(THETA0,AUG)
%      THETA0 - Initialization of Markov chain.  Each column
%               corresponds to a parameter and each row is a 
%               realization of the Markov chain.  If THETA0 is a 
%               row vector then it is considered the first element
%               of the Markov chain.  If THETA0 is a matrix then 
%               THETA0(end,:) is used in continuing the Markov
%               chain.  Thus, MCMC may be called successively to 
%               continue evaluation until satisfactory convergence is met.
%               If an InDir is given in the CONFIG.OPTIONS.MCMC field, then
%               the input is ignored and the parameter guess is taken from 
%               the last line of the highest ranking mat-file in the InDir. 
%
%               Note: THETA0 may be an array of structs.  In this case the 
%               fields of Theta0 are treated as "parameters" but need not 
%               be numeric.  If THETA0(end) is taken as the initial parameter
%               struct. 
%                
% 
%         AUG - Variable to store auxilliary information passed to and from 
%               likelihood and prior functions.     
% 
%
%	
%See also Configure.m
%

global CONFIG;
global LOGS;

global ROLLBACK;
ROLLBACK = false;

persistent TRACE_FIG;

%Determine whether messages will be stored into LOGS with HandleMessages

// if (nargin < 1)
//     if(HM), HandleMessages(3, 'At least 1 input argument is required.'); end;
// 	error('At least 1 input argument is required.')
// elseif (nargin > 2)
//     if(HM), HandleMessages(3, 'Too many input arguments.'); end;
// 	error('Too many input argmuments.');
// else
// 	if (nargin == 1)
// 		NEW_AUG = [];
// 		NEW_LOGBF = [];
// 	elseif (nargin == 2)	
// 		NEW_AUG = varargin{1};
// 		NEW_LOGBF = [];
// 	end
// end



TRACE_IGNORE = [];
PSTRUCT = false;



// if ~isempty(CONFIG.OPTIONS.MCMC.InDir) %Read parameter guesses from file -- ignore guesses given in input.
// 	if ~exist(CONFIG.OPTIONS.MCMC.InDir,'dir')
//         if(HM), HandleMessages(3, 'InDir could not be opened.'); end;
// 		error('InDir could not be opened.');
//     else
// 		S = dir(fullfile(CONFIG.OPTIONS.MCMC.InDir,'out_*.mat'));
// 		if (isempty(S))
//             if(HM), HandleMessages(3, 'Can''t get guesses -- no OUT files in InDir.'); end;
// 			error('Can''t get guesses -- no OUT files in InDir.');
// 		end
// 		
// 		%Find mat file with largest OUTNUM.
// 		N = strvcat(S.name);
// 		N = strrep(strrep(cellstr(N),{'out_'},{''}),{'.mat'},{''});
// 		N = str2num(strvcat(N));
//         
// 		[~,k] = max(N);
// 		FILE = fullfile(CONFIG.OPTIONS.MCMC.InDir,S(k).name);
// 		disp(['Loading last guesses from ',FILE]);
// 		Temp = load(FILE,'X');
//         
//         Xinit = Temp.X;
// 		NEW_AUG = Temp.AUG;
// 		NEW_LOGBF = Temp.LOGBF;
// 		%Xinit = Xinit(end,:);
// 	end
// end


// if isstruct(Xinit)
// 	PSTRUCT = true;
// 	FIELDNAMES = fieldnames(Xinit);
// 	LEGENDNAMES = strrep(FIELDNAMES,'_','\_');
// 	Xinit = [Xinit(:)]; %Ensure column depth orientation...	
// 	if strcmp(CONFIG.OPTIONS.MCMC.Trace,'on')
// 		temp = struct2cell(Xinit);
// 		TRACE_IGNORE = false(1,length(temp));
// 		
// 		for i = 1:length(temp)
// 			if ~isnumeric(temp{i})
//                 if(HM), HandleMessages(4, 'Trace can not be constructed with non-numeric parameters.'); end;
// 				warning('Trace can not be constructed with non-numeric parameters.')
// 				TRACE_IGNORE(i) = 1;
// 				disp(' ');
// 				disp(['Ignoring ''',FIELDNAMES{i},''' in trace.']);
// 				disp(' ');
// 			elseif ~isequal(size(temp{i}),[1 1])
// 				if(HM), HandleMessages(4, 'Trace can not be constructed with vector-valued parameters.'); end;
//                 warning('Trace can not be constructed with vector-valued parameters.')
// 				TRACE_IGNORE(i) = 1;
// 				disp(' ');
// 				disp(['Ignoring ''',FIELDNAMES{i},''' in trace.']);
// 				disp(' ');
// 			end
// 		end
// 		
// 		if ~any(TRACE_IGNORE) %Just to simplify checking later on we set to empty (defined off state)...
// 			TRACE_IGNORE = [];
// 		end
// 	end
// end


// %The default is to take one cycle with all parameters in the same block.
// if isempty(CONFIG.OPTIONS.MCMC.Blocks)
// 	if PSTRUCT
// 		FIELDNAMES = fieldnames(Xinit);
// 		CONFIG.OPTIONS.MCMC.Blocks = ones(1,length(FIELDNAMES));
// 	else
// 		CONFIG.OPTIONS.MCMC.Blocks = ones(1,size(Xinit,2));
// 	end
// end
// 
// if all(all(CONFIG.OPTIONS.MCMC.Blocks == 1))
// 	BLOCKING = false;
// else
// 	BLOCKING = true;
// end


%Take MAXT to be the number of complete block cycles.
if (size(Xinit,1) > 1)
	MAXT = size(Xinit,1) + CONFIG.OPTIONS.MCMC.MaxT*size(CONFIG.OPTIONS.MCMC.Blocks,1);
else	
	MAXT = CONFIG.OPTIONS.MCMC.MaxT*size(CONFIG.OPTIONS.MCMC.Blocks,1);
end


// if PSTRUCT
// 	X = repmat(Xinit(end),CONFIG.OPTIONS.MCMC.WriteOut,size(Xinit,2)); %Preallocate the X.
// else
// 	X = nan(CONFIG.OPTIONS.MCMC.WriteOut,size(Xinit,2)); %Preallocate the X.
// end

AUG = [];
LOGBF = [];

// if ~isempty(NEW_AUG)
// 	AUG = nan(CONFIG.OPTIONS.MCMC.WriteOut,max(1,size(NEW_AUG,2))); %Note that we pre-allocate to the current size in WriteOut
// 							%This may be larger or smaller than the past WriteOut size.
// end

if ~isempty(NEW_LOGBF)
	LOGBF = nan(CONFIG.OPTIONS.MCMC.WriteOut,max(1,size(NEW_AUG,2)));
end


t = size(Xinit,1); %Time (or iteration number).
row = t; %This is an index into X but will not always equal t.  For example row gets reset
	%to 1 after a WriteOut.

%X
%Xinit
if PSTRUCT
	X(1:length(Xinit)) = Xinit;
else

	X(1:size(Xinit,1),:) = Xinit; %Thus, the num. of rows of X is WriteOut or the size of Xinit.
			%whichever is bigger.
end

%If we loaded any values from a Mat-file, copy those values "new" values into the AUG and LOGBF... 
if ~isempty(NEW_AUG)
	AUG(1:size(Xinit,1),:) = NEW_AUG;
end
%LOGBF
if ~isempty(NEW_LOGBF)
	size(LOGBF)
	size(NEW_LOGBF)
	size(Xinit,1)
	LOGBF(1:size(Xinit,1),:) = NEW_LOGBF;
end

if isempty(LOGBF)
	LOGBF = nan(size(X,1),1);
	NEW_LOGBF = LOGBF;
end

  
// %Xinit = Xinit(:);
// %Xinit = Xinit';
// OUTFLAG = 'start';
OUTPUT_DATA([],[],[],OUTFLAG); %Initialize the data directory and files -- set OUTNUM persistent variable.

ACCEPT = 1;

if CONFIG.OPTIONS.MCMC.Display
	DISP_FINAL = true;
else
	DISP_FINAL = false;
end



ACCEPT_COUNT = 0;
TRY_COUNT = 0;
BLIKE = [];
BPRIOR = [];
BPROP = [];
ALIKE = [];
APRIOR = [];
APROP= [];

while true

	if (strcmp(OUTFLAG,'start')) %Initialize.
		if CONFIG.OPTIONS.MCMC.Display
			fprintf(1,'Iter\t');
			fprintf(1,'Block\t');
			fprintf(1,'Alpha\t');
			fprintf(1,'Result\n');
			fprintf(1,'%d\t',t);
			fprintf(1,'%d\t',1);
			fprintf(1,'%g\t',1);
			fprintf(1,'Initial\n');
		end
		
		if CONFIG.OPTIONS.MCMC.Trace %Plot a trace...

			if isempty(TRACE_FIG)
				TRACE_FIG = figure;
				set(TRACE_FIG,'Tag','MCMCTOOLS_TRACE');
            end
            
            if PSTRUCT

                if ~isempty(TRACE_IGNORE)
                    TEMP = struct2cell(X(1:row));
                    TEMP = TEMP(~TRACE_IGNORE,:);
                    TEMP = cell2mat(TEMP)';
                    the_line = line(t-row+1:t,TEMP);
                   

                    legend(LEGENDNAMES{~TRACE_IGNORE},-1);
                else
                    TEMP = struct2cell(X(1:row));
                    TEMP = cell2mat(TEMP)';
                    the_line = line(t-row+1:t,TEMP);
                    legend(LEGENDNAMES{:},-1);
                end

            else
                the_line = line(t-row+1:t,X(1:row,:));
            end

			drawnow;
		end
		
		OUTFLAG = 'in';
		
	else %Inside an iteration...	
		
		for BLOCK = 1:size(CONFIG.OPTIONS.MCMC.Blocks,1)
		
		
			BB = CONFIG.OPTIONS.MCMC.Blocks(BLOCK,:);
			
			%Get indices of parameters which should be changing.
			pset = (BB ~= 0);
			C = BB(pset); %The number of cycles for this parameter set.
			
			if (isempty(C))
                if(HM), HandleMessages(3, 'Parameter block cannot be empty. Check your MCMC block settings.'); end;
				error('Parameter block can not be empty. Check your MCMC block settings.');
			end
			C = C(1);
			
			
			
			for CYCLE = 1:C %No need to index here.
		
				
				%Get a random occurence from the proposal distribution.
				%
				%Note that in Gibbs sampling the proposal should be the 
				%posterior distribution (i.e. you should be sampling 
				%directly from the posterior distribution).  Gibbs
				%sampling may not seem useful at first.  The idea
				%is that for some very complicated posteriors that
				%we have an analytic expression for may still be 
				%difficult to calculate expectations over etc.
				%Gibbs sampling allows us to get means and variances
				%of these really complicated posteriors we don't know.
				%
				while (true)
					ROLLBACK = false;
					for rolls=1:1
						Y = feval(CONFIG.FUNCTIONS.Proposal,'rand',X(row,:));

						NEW_AUG = [];


						if ~isempty(AUG)
							NEW_AUG = AUG(row,:);
						end



						TRY_COUNT = TRY_COUNT+1;



						if PSTRUCT
							if ~isstruct(Y)
                                if(HM), HandleMessages(3, 'Proposal must be a struct if parameter struct is used.'); end;
								error('Proposal must be a struct if parameter struct is used.')
							end
						end	

						switch CONFIG.OPTIONS.MCMC.Method
                            case 'mh' %Metropolis-Hastings

                                %AUG is a variable containing auxilliary information.
                                %The original reason for its intorduction was
                                %to keep track of data augmentation.
                                %

                                if(ACCEPT ==1 &&~isempty(ALIKE)&&~isempty(APRIOR))
                                                %disp('Setting accepted values to current values.');
                                                BLIKE = ALIKE;
                                                BPRIOR = APRIOR;
                                            end

                                            %disp('Evaluating Likelihood.');
                                [ALIKE,NEW_AUG] = feval(CONFIG.FUNCTIONS.Likelihood,Y,NEW_AUG);

                                if (ROLLBACK)
                                    if(HM)
                                        HandleMessages(2, 'Rollback!');
                                    else
                                        disp('Rollback!');
                                    end
                                    break;%Breaks out of inner rollback loop
                                end
                                            %disp('Evaluating Prior distribution');
                                [APRIOR,NEW_AUG] = feval(CONFIG.FUNCTIONS.Prior,Y,NEW_AUG);

                                if (ROLLBACK)
                                    if(HM)
                                        HandleMessages(2, 'Rollback!');
                                    else
                                        disp('Rollback!');
                                    end
                                    break;%Breaks out of inner rollback loop
                                end


                                if (~CONFIG.OPTIONS.MCMC.AssertSymmetricProposal)
                                    APROP = feval(CONFIG.FUNCTIONS.Proposal,'pdf',X(row,:),Y);
                                end


                                if (ROLLBACK)
                                    if(HM)
                                        HandleMessages(2, 'Rollback!');
                                    else
                                        disp('Rollback!');
                                    end
                                    break;%Breaks out of inner rollback loop
                                end

                                if any(isnan(ALIKE))
                                    if(HM), HandleMessages(3, ['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Likelihood,'.']); end;
                                    error(['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Likelihood,'.']);
                                end

                                if any(isnan(APRIOR))
                                    if(HM), HandleMessages(3, ['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Prior,'.']); end;
                                    error(['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Prior,'.']);
                                end

                                if any(isnan(APROP))
                                    if(HM), HandleMessages(3, ['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Proposal,'.']); end;
                                    error(['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Proposal,'.']);
                                end

                                if (isempty(BLIKE) || CONFIG.OPTIONS.MCMC.ReCalc)
                                        %disp('BLIKE is empty, filling it in...');
                                        [BLIKE,NEW_AUG] = feval(CONFIG.FUNCTIONS.Likelihood,X(row,:),NEW_AUG);
                                end


                                if (ROLLBACK)
                                    if(HM)
                                        HandleMessages(2, 'Rollback!');
                                    else
                                        disp('Rollback!');
                                    end
                                    break;%Breaks out of inner rollback loop
                                end

                                if (isempty(BPRIOR) || CONFIG.OPTIONS.MCMC.ReCalc)
                                    %disp('BPRIOR is empty, filling it in...');
                                    [BPRIOR,NEW_AUG] = feval(CONFIG.FUNCTIONS.Prior,X(row,:),NEW_AUG);
                                end

                                if (ROLLBACK)
                                    if(HM)
                                        HandleMessages(2, 'Rollback!');
                                    else
                                        disp('Rollback!');
                                    end
                                    break;%Breaks out of inner rollback loop
                                end

                                if any(isnan(BLIKE))
                                    if(HM), HandleMessages(3, ['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Likelihood,'.']); end;
                                    error(['Nan''s detected in evaluation of ', CONFIG.FUNCTIONS.Likelihood,'.']);
                                end

                                if any(isnan(BPRIOR))
                                    if(HM), HandleMessages(3, ['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Prior,'.']); end;
                                    error(['Nan''s detected in evaluation of ',Prior,'.']);
                                end

                                if (~CONFIG.OPTIONS.MCMC.AssertSymmetricProposal) %Proposal is not symmetric.
                                    BPROP = feval(CONFIG.FUNCTIONS.Proposal,'pdf',Y,X(row,:));
                                end

                                if (ROLLBACK)
                                    if(HM)
                                        HandleMessages(2, 'Rollback!');
                                    else
                                        disp('Rollback!');
                                    end
                                    break;%Breaks out of inner rollback loop
                                end

                                if isnan(BPROP)
                                    if(HM), HandleMessages(3, ['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Proposal,'.']); end;
                                    error(['Nan''s detected in evaluation of ',CONFIG.FUNCTIONS.Proposal,'.']);
                                end
                                lastwarn('');

                                if (any(size(ALIKE)~=size(BLIKE)))
                                    s1 = sprintf('%dx',size(ALIKE));
                                    s1 = s1(1:end-1);
                                    s2 = sprintf('%dx',size(BLIKE));
                                    s2 = s2(1:end-1);

                                    if(HM)
                                        HandleMessages(2, ['ALIKE [',s1,'] and BLIKE [',s2,']  not same size. Using what is available.']);
                                    else
                                        disp(['ALIKE [',s1,'] and BLIKE [',s2,']  not same size. Using what is available.']);
                                    end

                                    mn = min(size(ALIKE,1),size(BLIKE,1));
                                    TA = log(ALIKE(1:mn,:)./BLIKE(1:mn,:));
                                        else
                                                if(CONFIG.OPTIONS.MCMC.UseLogLike)
                                                    TA = ALIKE-BLIKE;
                                                else
                                                    TA = log(ALIKE./BLIKE);
                                                end
                                end

                                if ~isempty(lastwarn)
                                    if(HM)
                                        HandleMessages(2, ['Zero likelihood detected in evaluation of ', CONFIG.FUNCTIONS.Likelihood,'.']);
                                    else
                                        disp(['Zero likelihood detected in evaluation of ', CONFIG.FUNCTIONS.Likelihood,'.']);
                                    end
                                end

                                lastwarn('');

                                if (~CONFIG.OPTIONS.MCMC.AssertSymmetricProposal)
                                    TB = log(APROP./BPROP);
                                else
                                    TB = 0; %Assume factor APROP/BPROP simplified to 1 due to symmetry.
                                end

                                if ~isempty(lastwarn)
                                    if(HM)
                                        HandleMessages(2, ['Zero likelihood detected in evaluation of ', CONFIG.FUNCTIONS.Proposal,'.']);
                                    else
                                        disp(['Zero likelihood detected in evaluation of ', CONFIG.FUNCTIONS.Proposal,'.']);
                                    end
                                end

                                lastwarn('');
                                TC = log(APRIOR./BPRIOR);
                                if ~isempty(lastwarn)
                                    if(HM)
                                        HandleMessages(2, ['Zero likelihood detected in evaluation of ', CONFIG.FUNCTIONS.Prior,'.']);
                                    else
                                        disp(['Zero likelihood detected in evaluation of ', CONFIG.FUNCTIONS.Prior,'.']);
                                    end
                                end

                                NEW_LOGBF = sum([TA(:);TB(:);TC(:)]); 
                                alpha =  min(1,exp(NEW_LOGBF));
                                        %disp('End of an iteration.');
                                        %X(row)



                            case 'm' %Metropolis - use when full posterior can be evaluated 
                                %Note: usually the sticking part is in the denominator
                                %in Bayes rule (i.e. the P(B) in P(A|B) = P(B|A)*P(A)/P(B))
                                %If we don't have P(B) then we use Metrop. Hastings.  If we
                                %do have this then we can use Metropolis.
                                %

                                A = feval(CONFIG.FUNCTIONS.Posterior,Y);

                                if isnan(A)
                                    if(HM), HandleMessages(3, 'Nan''s detected'); end;
                                    error('Nan''s detected');
                                end


                                B = feval(CONFIG.OPTIONS.Posterior,X(row,:));	

                                if isnan(B)
                                    if(HM), HandleMessages(3, 'Nan''s detected'); end;
                                    error('Nan''s detected');
                                end

                                alpha = min(1,A./B);

                            case 'gibbs' %Gibbs -- use when sampling directly from the posterior.

                                alpha = 1;  %Candidate is always accepted!
                            otherwise
                                if(HM), HandleMessages(3, ['Unknown method ' CONFIG.OPTIONS.MCMC.Method '.']); end;
                                error(['Unknown method ' CONFIG.OPTIONS.MCMC.Method '.'])
						end
					end %END ROLLBACK LOOP
					
					if (~ROLLBACK)
						break; %Breaks out of enclosing WHILE(TRUE) loop
					end	
				
				end %END OF WHILE(TRUE) LOOP for ROLLBACK
				
				if (alpha > rand(1))
					ACCEPT_COUNT = ACCEPT_COUNT+1;
					if CONFIG.OPTIONS.MCMC.Display
						fprintf(1,'%d\t',t);
						fprintf(1,'%d\t',BLOCK);
						fprintf(1,'%g\t',alpha);
						fprintf(1,'Accept\n');
					end
					
					ACCEPT = 1;
					
					if PSTRUCT
						if BLOCKING
							temp = struct2cell(Y);
							temp2 = struct2cell(X(row));
							[temp{~pset}] = deal(temp2{~pset});
							X(row+1) = cell2struct(temp,FIELDNAMES);
						else
							X(row+1) = Y;
						end
					else
						X(row+1,pset) = Y(1,pset);
						X(row+1,~pset) = X(row,~pset);
					end
					
					if ~isempty(NEW_AUG)
						AUG(row+1,:) = NEW_AUG(:)';
					end
					
					
					LOGBF(row+1,:) = NEW_LOGBF(:)';
					
					if (CONFIG.OPTIONS.MCMC.kdetrack)
						%Copy "latest" kde queue directory to a subdirecotry withname having the row as a tracking number.  
						SOURCE = CONFIG.OPTIONS.MCMC.kdequeue;
						
						if (exist(SOURCE,'dir'))
							DESTINATION = fullfile(CONFIG.OPTIONS.MCMC.kdesave,['output_',num2str(row+1)]);
							if (~exist(CONFIG.OPTIONS.MCMC.kdesave,'dir'))
								mkdir(CONFIG.OPTIONS.MCMC.kdesave);
							end

							if (~exist(DESTINATION,'dir'))
								mkdir(DESTINATION);
							end
							
							movefile(SOURCE,DESTINATION);
						end
					end
				
				else
					if CONFIG.OPTIONS.MCMC.Display
						fprintf(1,'%d\t',t);
						fprintf(1,'%d\t',BLOCK);
						fprintf(1,'%g\t',alpha);
						fprintf(1,'Reject\n');
					end
					
					ACCEPT = 0;
							
					X(row+1,:) = X(row,:);
					if (isnan(LOGBF(row,:)) && (row == 1))
						%We make a special exception for the case
						%when the first parameter guess results in a
						%rejection.  In that case we set the initial 
						%log bayes to the log bayes of the rejected guess.
						%
						LOGBF(row+1,:) = NEW_LOGBF(:)';
					else
						LOGBF(row+1,:) = LOGBF(row,:);
					end
					
					if ~isempty(AUG)
						AUG(row+1,:) = AUG(row,:);
					end
					
					if (CONFIG.OPTIONS.MCMC.kdetrack)
						%Due to rejection, we ...
						%Copy "previous" kde queue directory to a subdirecotry withname having the row as a tracking number.  
						SOURCE = fullfile(CONFIG.OPTIONS.MCMC.kdesave,['output_',num2str(row)]);
						if (exist(SOURCE,'dir'))
							DESTINATION = fullfile(CONFIG.OPTIONS.MCMC.kdesave,['output_',num2str(row+1)]);

							if (~exist(CONFIG.OPTIONS.MCMC.kdesave,'dir'))
								mkdir(CONFIG.OPTIONS.MCMC.kdesave);
							end
							%if (~exist(DESTINATION,'dir'))
	%							mkdir(DESTINATION);
	%						end
							movefile(CONFIG.OPTIONS.MCMC.kdequeue,DESTINATION);
						end
					end
				end
				
				
				if CONFIG.OPTIONS.MCMC.Trace
					%TRACE_FIG = findobj('Type','figure','Tag','MCMCTOOLS_TRACE');
					
					if isempty(TRACE_FIG)
						TRACE_FIG = figure;
						%disp(['CREATING TRACE_FIG=',num2str(TRACE_FIG)]);
						
						set(TRACE_FIG,'Tag','MCMCTOOLS_TRACE');
					else
						%disp(['CURRENTLY TRACE_FIG=',num2str(TRACE_FIG)]);
						TRACE_FIG = TRACE_FIG(1); %Get the first if multiple figures have same tag.
                    end
					
					
					if ~exist('the_line')  %This wouold happen if the line got deleted somehow...
						figure(TRACE_FIG);
						%Trace the initial condition vector
						%or matrix.
						if PSTRUCT
							if ~isempty(TRACE_IGNORE)
							
								TEMP = struct2cell(X(1:row));
								TEMP = TEMP(~TRACE_IGNORE,:);
								TEMP = cell2mat(TEMP)';
								figure(TRACE_FIG);
								the_line = line(t-row+1:t,TEMP);
								%legend(FIELDNAMES{~TRACE_IGNORE},-1);
							else
								TEMP = struct2cell(X(1:row));
								TEMP = cell2mat(TEMP)';
								figure(TRACE_FIG);
								the_line = line(t-row+1:t,TEMP);
								%legend(FIELDNAMES{:},-1);
							end
							
						else
							the_line = line(t-row+1:t,X(1:row,:));
						end


					else
						figure(TRACE_FIG);
						%disp(['TRACE_FIG=',num2str(TRACE_FIG)]);
						
						if PSTRUCT
							
							if ~isempty(TRACE_IGNORE)
								TEMP = struct2cell(X(row));
								TEMP = TEMP(~TRACE_IGNORE,:);
								TEMP = cell2mat(TEMP)';
							else
								
								TEMP = struct2cell(X(row));
								TEMP = cell2mat(TEMP)';
								
							end
						else
							
							TEMP = X(row,:);
							
						end
						
						for i = 1:length(the_line)
							if PSTRUCT
	
								old_x = get(the_line(i),'XData');
								new_x = [old_x,t];
								old_y = get(the_line(i),'YData');
								new_y = [old_y,TEMP(i)'];
								set(the_line(i),'XData',new_x,'YData',new_y);

							else
								old_x = get(the_line(i),'XData');
								new_x = [old_x,t];
								old_y = get(the_line(i),'YData');
								new_y = [old_y,TEMP(i)'];
								set(the_line(i),'XData',new_x,'YData',new_y);
							end
						end


					end
					drawnow;
					
				end
				
				row = row+1;
				t = t+1; %We have just completed 1 more iteration.
				
				
				
				if (t >= MAXT)
					
					if CONFIG.OPTIONS.MCMC.Display
                        if(HM)
                            HandleMessages(2, 'Max iterations reached.')
                        else
                            disp('Max iterations reached.');
                        end
					end
					%Cut off the trailing entries (NAN if numeric and fields if PSTRUCT)
					% if any.
					X = X(1:row,:);
					LOGBF = LOGBF(1:row,:);
					if ~isempty(AUG)
						AUG = AUG(1:row,:);
						
						OUTPUT_DATA(X,AUG,LOGBF,OUTFLAG);
					else
						OUTPUT_DATA(X,[],LOGBF,OUTFLAG);
					end
					
					RATIO = ACCEPT_COUNT/TRY_COUNT;
					if DISP_FINAL
                        if(HM)
                            HandleMessages(2, ['Acceptance ratio: ' num2str(RATIO)]);
                        else
                            fprintf('Acceptance ratio: %f',RATIO);		
                        end
					end

					return;
				else %Continue ...
                    %Iterate the logs
					HandleMessages('iterate');
                    
					if ((mod(t-1,CONFIG.OPTIONS.MCMC.WriteOut) == 0) && ~isempty(CONFIG.OPTIONS.MCMC.OutDir))
						%disp('Outputting parameter guesses.');
						
						if ~isempty(AUG)
							OUTPUT_DATA(X(1:row-1,:),AUG(1:row-1,:),LOGBF(1:row-1,:),OUTFLAG);
						else
							
							OUTPUT_DATA(X(1:row-1,:),[],LOGBF(1:row-1,:),OUTFLAG);
						end
						
						if PSTRUCT
							
							temp = X(row-1,1);
							F = fieldnames(temp);
							V = repmat({nan},length(F),1);
							
							X = repmat(cell2struct(V,F),size(X,1),1);
							
							X(1) = temp;
						else
							temp = X(row-1,:);
							X(:) = nan;
							X(1,:) = temp;
						end
						
						if ~isempty(AUG)
							temp = AUG(row-1,:);
							AUG(:) = nan;
							AUG(1,:) = temp;
						else
							AUG = [];
						end
						
						row = 1;
					end
					
					
					
				end
	
			end
		end
		
		
		
	end
	
end



function OUTPUT_DATA(X,AUG,LOGBF,varargin)
%OUTPUT_DATA append (possibly) parameter guesses to the output file.
%  OUTPUT_DATA is responsible for implementing the output options
%  using only the input FLAG.
global CONFIG;
global LOGS;

persistent OUTNUM;

if(CONFIG.OPTIONS.DEBUG.HandleMessages)
    HM = true;
else
    HM = false;
end

if isempty(CONFIG.OPTIONS.MCMC.OutDir) %Don't do anything.
    
    return;		
else
	FLAG = '';
	if (nargin == 4)
		FLAG = varargin{1};
    end
	
	
	switch FLAG
	
	case 'start'
		[A,B] = fileparts(CONFIG.OPTIONS.MCMC.OutDir);
		
		if isempty(A)
			A = pwd;
		end
		
		%Check to see if DeleteOnStart is 'yes'.
		if (CONFIG.OPTIONS.MCMC.DeleteOnStart) %Delete output files currently in the
			%directory.
			
			if exist(fullfile(A,B),'dir')
                if(HM)
                    HandleMessages(2, 'Deleting .mat files in OutDir...');
                else
                    disp('Deleting .mat files in OutDir...');
                end
				delete(fullfile(CONFIG.OPTIONS.MCMC.OutDir,'*_out_*.mat'));
                if(CONFIG.OPTIONS.GRAPHICS.DeleteGraphicsOnStart)
                    if(HM)
                        HandleMessages(2, 'Deleting graphics files in OutDir...');
                    else
                        disp('Deleting graphics files in OutDir...');
                    end
                    delete(fullfile(CONFIG.OPTIONS.MCMC.OutDir, [CONFIG.META.mcmcrunnickname '_*.' CONFIG.OPTIONS.GRAPHICS.graphicsfmt]));
                end
			else %Directory does not exist
                if(HM)
                    HandleMessages(2, 'Creating OutDir...');
                else
                    disp('Creating OutDir...');
                end
				[Success,MSG] = mkdir(A,B);
				if ~Success
					warning(MSG);
				end
			end
			
			OUTNUM = 1;
		else		
			if CONFIG.OPTIONS.MCMC.Append
				%Determine largest output number of matfiles in this dir.
				if ~exist(fullfile(A,B),'dir')
					[Success,MSG] = mkdir(A,B); %Create the directory.
					if ~Success
						warning(MSG);
					end
				else %Directory exists.
					S = dir(fullfile(CONFIG.OPTIONS.MCMC.InDir,[CONFIG.META.mcmcrunnickname '_out_*.mat']));
					
					N = strvcat(S.name);
					

					if isempty(N)
						OUTNUM = 1;
					else
					
						N = strrep(strrep(cellstr(N),{'out_'},{''}),{'.mat'},{''});
						N = str2num(strvcat(N));
						
						OUTNUM = max(N)+1;
					end
				end
				
			else
				OUTNUM = 1; %We will just be overwriting the file out_1.mat.
			end
		end
		
		
		
		return; %Just return at this point!!
		
	case 'in'
        if(HM)
            HandleMessages(2, 'Writing out');
        else
            disp('Writing out');
        end
        
		mat_file = [CONFIG.META.mcmcrunnickname '_out_',num2str(OUTNUM),'.mat'];
		mat_file = fullfile(CONFIG.OPTIONS.MCMC.OutDir,mat_file);
		
		save(mat_file,'X','AUG','LOGBF', 'LOGS'); %Save data to a binary matfile within the OutDir called <mcmcrunname>_out_<OUTNUM>
		
		
		if CONFIG.OPTIONS.MCMC.Append
		        OUTNUM = OUTNUM+1;
		else %Overwrite...
			OUTNUM = 1;
		end
	end
	
	
	
end


