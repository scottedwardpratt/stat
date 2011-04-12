function mcmc_sens(PRESULT,PROPFUN,PROPPARAMS,LIKEFUN,LIKEPARAMS,PRIORFUN,PRIORPARAMS,[],OPTS);
%Measure sensitivity of the posterior distribution to parameter perturbations.
%
%

if (nargin < 7)
	error('Atleast 7 input arguments are required.')
elseif (nargin > 9)
	error('Too many input argmuments.');
else
	if (nargin == 7)
		NEW_AUG = [];
		OPTS = mcmcset;
	elseif (nargin == 8)	
		NEW_AUG = varargin{1};
		OPTS = mcmcset;
		
	elseif (nargin == 9)
		NEW_AUG = varargin{1};
		
		if ~isstruct(varargin{2})
			error('Options must be struct.');
		else
			OPTS = mcmcset(varargin{2});
		end
	end
end

TRACE_IGNORE = [];
PSTRUCT = logical(0);
if isstruct(Xinit)
	PSTRUCT = logical(1);
	FIELDNAMES = fieldnames(Xinit);
	Xinit = [Xinit(:)]; %Ensure column depth orientation...	
	if strcmp(OPTS.Trace,'on')
		temp = struct2cell(Xinit);
		TRACE_IGNORE = logical(zeros(1,length(temp)));
		
		for i = 1:length(temp)
			if ~isnumeric(temp{i})
				warning('Trace can not be constructed with non-numeric parameters.')
				TRACE_IGNORE(i) = 1;
				disp(' ');
				disp(['Ignoring ''',FIELDNAMES{i},''' in trace.']);
				disp(' ');
			elseif ~isequal(size(temp{i}),[1 1])
				warning('Trace can not be constructed with vector-valued parameters.')
				TRACE_IGNORE(i) = 1;
				disp(' ');
				disp(['Ignoring ''',FIELDNAMES{i},''' in trace.']);
				disp(' ');
			end
		end
		
		if ~any(TRACE_IGNORE) %Just to simplify checking later on...
			TRACE_IGNORE = [];
		end
	end
end


%The default is to take one cycle with all parameters in the same block.
if isempty(OPTS.Blocks)
	if PSTRUCT
		OPTS.Blocks = ones(1,length(FIELDNAMES));
	else
		OPTS.Blocks = ones(1,size(Xinit,2));
	end
end

if all(all(OPTS.Blocks == 1))
	BLOCKING = logical(0);
else
	BLOCKING = logical(1);
end





if ~isempty(OPTS.InDir) %Read parameter guesses from file -- ignore guesses given in input.
	if ~exist(OPTS.InDir,'dir')
		error('InDir could not be opened.');
	else
		
		S = dir(fullfile(OPTS.InDir,'out_*.mat'));
		if (length(S) == 0)
			error('Can''t get guesses -- no OUT files in InDir.');
		end
		
		%Find mat file with largest OUTNUM.
		N = strvcat(S.name);
		ind_out = find(strcmp(cellstr(N(:,1:4)),'out_'));
		N = str2num(N(ind_out,5:end-4));
		[junk,ind] = max(N);
		Xinit = getfield(load(fullfile(OPTS.InDir,S(ind_out(ind)).name),'X'),'X');
		NEW_AUG = getfield(load(fullfile(OPTS.InDir,S(ind_out(ind)).name),'AUG'),'AUG');
		
		
		%Xinit = Xinit(end,:);
	end
end



%Take MAXT to be the number of complete block cycles.
if (size(Xinit,1) > 1)
	
	MAXT = size(Xinit,1) + OPTS.MaxT*size(OPTS.Blocks,1);
	
else	
	MAXT = OPTS.MaxT*size(OPTS.Blocks,1);
end

if PSTRUCT
	X = repmat(Xinit(end),OPTS.WriteOut,size(Xinit,2)); %Preallocate the X.
else
	X = repmat(nan,OPTS.WriteOut,size(Xinit,2)); %Preallocate the X.
end

AUG = [];
if ~isempty(NEW_AUG)
	AUG = repmat(nan,OPTS.WriteOut,size(NEW_AUG,2));
end

t = size(Xinit,1); %Time (or iteration number).
row = t; %This is an index into X but will not always equal t.  For example row gets reset
	%to 1 after a WriteOut.
X(1:size(Xinit,1),:) = Xinit; %Thus, the num. of rows of X is OPTS.WriteOut or the size of Xinit.
			%whichever is bigger.


if ~isempty(NEW_AUG)
	AUG(1:size(Xinit,1),:) = NEW_AUG;
end

%Xinit = Xinit(:);
%Xinit = Xinit';
OUTFLAG = 'start';
OUTPUT_DATA([],[],OPTS,OUTFLAG); %Initialize the data directory and files -- set OUTNUM persistent variable.

ACCEPT = 1;

if ~strcmp(lower(OPTS.Display),'off')
	DISP_FINAL = logical(1);
else
	DISP_FINAL = logical(0);
end

ACCEPT_COUNT = 0;
TRY_COUNT = 0;

while 1	
	
	if (strcmp(OUTFLAG,'start')) %Initialize.
		if strcmp(OPTS.Display,'on')
			fprintf(1,'Iter\t');
			fprintf(1,'Block\t');
			fprintf(1,'Alpha\t');
			fprintf(1,'Result\n');
			fprintf(1,'%d\t',t);
			fprintf(1,'%d\t',1);
			fprintf(1,'%g\t',1);
			fprintf(1,'Initial\n');
		end
		
		if strcmp(OPTS.Trace,'on') %Plot a trace...
			%h1 = findobj('Type','figure','Tag','MCMCTOOLS_TRACE');

			% if isempty(h1)
% 				h1 = figure;
% 				set(h1,'Tag','MCMCTOOLS_TRACE');
% 			else
% 				h1 = h1(1); %Get the first if multiple figures have same tag.
% 				newplot;
% 			end
% 			

			newplot;
			h1 = gcf;
			set(h1,'Tag','MCMCTOOLS_TRACE');
			
			%if ~exist('the_line')
				%figure(h1);
				%Trace the initial condition vector
				%or matrix.
				if PSTRUCT
					
					if ~isempty(TRACE_IGNORE)
						TEMP = struct2cell(X(1:row));
						TEMP = TEMP(~TRACE_IGNORE,:);
						TEMP = cell2mat(TEMP)';
						the_line = line([t-row+1:t],TEMP);
						%text(t,TEMP(end,:)
						
						legend(FIELDNAMES{~TRACE_IGNORE},-1);
					else
						%TEMP = cell2mat(struct2cell(X(1:row)));
						TEMP = struct2cell(X(1:row));
						TEMP = cell2mat(TEMP)';
						the_line = line([t-row+1:t],TEMP);
						legend(FIELDNAMES{:},-1);
					end
					
				else
					the_line = line([t-row+1:t],X(1:row,:));
				end
				
				
			% else
% 				figure(h1);
% 				
% 				if PSTRUCT
% 					if ~isempty(TRACE_IGNORE)
% 						TEMP = struct2cell(X(1:row));
% 						TEMP = TEMP(~TRACE_IGNORE,:);
% 						TEMP = cell2mat(TEMP)';
% 						
% 					else
% 						TEMP = struct2cell(X(1:row));
% 						TEMP = cell2mat(TEMP)';
% 						
% 					end
% 				end
% 				
% 				for i = 1:length(the_line)
% 					if PSTRUCT
% 						old_x = get(the_line(i),'XData');
% 						new_x = [old_x,t];
% 						old_y = get(the_line(i),'YData');
% 						new_y = [old_y,TEMP(i)'];
% 						set(the_line(i),'XData',new_x,'YData',new_y);
% 						
% 					else
% 						old_x = get(the_line(i),'XData');
% 						new_x = [old_x,t];
% 						old_y = get(the_line(i),'YData');
% 						new_y = [old_y,TEMP(i)'];
% 						set(the_line(i),'XData',new_x,'YData',new_y);
% 					end
% 				end
% 				
% 				
% 			end
			drawnow;
		end
		
		OUTFLAG = 'in';
		
	else %Inside an iteration...	
		
	
		for BLOCK = 1:size(OPTS.Blocks,1)
		
		
			BB = OPTS.Blocks(BLOCK,:);
		
			%Get indices of parameters which should be changing.
			pset = (BB ~= 0);
			C = BB(pset); %The number of cycles for this parameter set.
			C = C(1);
			
			
			for CYCLE = 1:C %No need to index here.
		
				
				%Get a random occurence from the proposal distribution.
				Y = feval(proposal,'rand',X(row,:),[],propparams{:});
				
				NEW_AUG = [];
				if ~isempty(AUG)
					NEW_AUG = AUG(row,:);
				end
				
				TRY_COUNT = TRY_COUNT+1;
				
				
				
				if PSTRUCT
					if ~isstruct(Y)
						error('Proposal must be a struct if parameter struct is used.')
					end
				end	
				
				switch OPTS.Method
				case 'mh' %Metropolis-Hastings
					
					%A = [feval(posterior,Y,postparams{:})*feval(proposal,'pdf',X(b,:),Y,propparams{:})];
					
					%AUG is a variable containing auxilliary information.
					%The original reason for its intorduction was
					%to keep track of data augmentation.
					%
					
					[ALIKE,NEW_AUG] = feval(likelihood,Y,NEW_AUG,likeparams{:});
					
					[APRIOR,NEW_AUG] = feval(prior,Y,NEW_AUG,priorparams{:});
					APROP = feval(proposal,'pdf',X(row,:),Y,propparams{:});
					
					if any(isnan(ALIKE))
						error(['Nan''s detected in evaluation of ',likelihood,'.']);
					end

					if any(isnan(APRIOR))
						error(['Nan''s detected in evaluation of ',prior,'.']);
					end

					if any(isnan(APROP))
						error(['Nan''s detected in evaluation of ',proposal,'.']);
					end

					
					
					if strcmp(OPTS.ReCalc,'on')
						%If the previous proposal value of the chain was accepted						
						%then we need to recalculate this.  Otherwise it is wasting
						%computational effort to recalculate these!!
						%
						%Note however that in cases such as data augmentation
						%this slows down mixing and although each iteration 
						%takes longer overall convergence is faster.
						%
						
						[BLIKE,NEW_AUG] = feval(likelihood,X(row,:),NEW_AUG,likeparams{:});
						[BPRIOR,NEW_AUG] = feval(prior,X(row,:),NEW_AUG,priorparams{:});
						% if strcmp(likelihood,'lfit_like')
% 							disp('Here')
% 							row
% 							X(row,:)
% 							%X(row+1,:)
% 							X
% 							
% 						end
						if any(isnan(BLIKE))
							
							error(['Nan''s detected in evaluation of ',likelihood,'.']);
						end

						if any(isnan(BPRIOR))
							error(['Nan''s detected in evaluation of ',prior,'.']);
						end
					end
					
					
					BPROP = feval(proposal,'pdf',Y,X(row,:),propparams{:});

					if isnan(BPROP)
						error(['Nan''s detected in evaluation of ',proposal,'.']);
					end
					lastwarn('');
					TA = log(ALIKE./BLIKE);
					if ~isempty(lastwarn)
						disp(['Zero likelihood detected in evaluation of ', likelihood,'.']);
					end
					
					lastwarn('');
					TB = log(APROP./BPROP);
					if ~isempty(lastwarn)
						disp(['Zero likelihood detected in evaluation of ', proposal,'.']);
					end
					
					lastwarn('');
					TC = log(APRIOR./BPRIOR);
					if ~isempty(lastwarn)
						disp(['Zero likelihood detected in evaluation of ', prior,'.']);
					end
					
					LOGBF = sum([TA(:);TB(:);TC(:)]); 
					
 					alpha =  min(1,exp(LOGBF));



				case 'm' %Metropolis
					A = feval(posterior,Y,postparams{:});

					if isnan(A)
						error('Nan''s detected');
					end


					B = feval(posterior,X(row,:),postparams{:});	

					if isnan(B)
						error('Nan''s detected');
					end

					alpha = min(1,A./B);

				case 'gibbs' %Gibbs

					alpha = 1;  %Candidate is always accepted!
				otherwise
					error('Unknown method.')
				end


				
				
				if (alpha > rand(1))
					ACCEPT_COUNT = ACCEPT_COUNT+1;
					if strcmp(OPTS.Display,'on')
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
					
				else
					if strcmp(OPTS.Display,'on')
						fprintf(1,'%d\t',t);
						fprintf(1,'%d\t',BLOCK);
						fprintf(1,'%g\t',alpha);
						fprintf(1,'Reject\n');
					end
					
					ACCEPT = 0;
							
					X(row+1,:) = X(row,:);
					
					if ~isempty(AUG)
						AUG(row+1,:) = AUG(row,:);
					end
				end
				
				
				if strcmp(OPTS.Trace,'on')
					h1 = findobj('Type','figure','Tag','MCMCTOOLS_TRACE');
					
					if isempty(h1)
						h1 = figure;
						set(h1,'Tag','MCMCTOOLS_TRACE');
					else
						h1 = h1(1); %Get the first if multiple figures have same tag.
					end
					
					% if ~exist('the_line')
% 						figure(h1);
% 						the_line = line([t-row+1:t],X(1:row,:));
% 						
% 					else
% 						figure(h1);
% 						
% 						for i = 1:length(the_line)
% 							old_x = get(the_line(i),'XData');
% 							new_x = [old_x,t];
% 							old_y = get(the_line(i),'YData');
% 							new_y = [old_y,X(row,i)'];
% 							set(the_line(i),'XData',new_x,'YData',new_y);
% 						end	
% 						
% 					end
% 					drawnow;
% 					
					
					
					if ~exist('the_line')  %This wouold happen if the line got deleted somehow...
						figure(h1);
						%Trace the initial condition vector
						%or matrix.
						if PSTRUCT
							if ~isempty(TRACE_IGNORE)
							
								TEMP = struct2cell(X(1:row));
								TEMP = TEMP(~TRACE_IGNORE,:);
								TEMP = cell2mat(TEMP)';
								the_line = line([t-row+1:t],TEMP);
								%legend(FIELDNAMES{~TRACE_IGNORE},-1);
							else
								TEMP = struct2cell(X(1:row));
								TEMP = cell2mat(TEMP)';
								
								the_line = line([t-row+1:t],TEMP);
								%legend(FIELDNAMES{:},-1);
							end
							
						else
							the_line = line([t-row+1:t],X(1:row,:));
						end


					else
						figure(h1);
						
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
					
					if strcmp(OPTS.Display,'on')
						disp('Max iterations reached.');
					end
					%Cut off the trailing entries (NAN if numeric and fields if PSTRUCT)
					% if any.
					X = X(1:row,:);
					
					
					if ~isempty(AUG)
						AUG = AUG(1:row,:);
						OUTPUT_DATA(X,AUG,OPTS,OUTFLAG);
					else
						OUTPUT_DATA(X,[],OPTS,OUTFLAG);
					end
					
					RATIO = ACCEPT_COUNT/TRY_COUNT;
					if DISP_FINAL
						disp(sprintf('Acceptance ratio: %f',RATIO));		
					end

					return;
				else %Continue ...
					if ((mod(t-1,OPTS.WriteOut) == 0) & ~isempty(OPTS.OutDir))
						%disp('Outputting parameter guesses.');
						
						if ~isempty(AUG)
							OUTPUT_DATA(X(1:row-1,:),AUG(1:row-1,:),OPTS,OUTFLAG);
						else
							OUTPUT_DATA(X(1:row-1,:),[],OPTS,OUTFLAG);
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



	       
