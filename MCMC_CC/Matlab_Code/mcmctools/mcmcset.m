function S = mcmcset(varargin)
%MCMCSET Options structure for MCMC.
%
%      Append - Add output file to OutDir during a call to MCMC
%       	[ {yes} no].  If set to 'no' a single output file
%       	will be overwritten at each WriteOut. 
%
%      Blocks - Matrix of integers indicating blocking as well as the number 
%       	of iterations taken for each parameter block.
%       	[ pos. integer matrix ]
%
%DeleteOnStart - If 'yes' the OutFile (if it exists) is overwritten
%               at the start of a call to MCMC. The output files  
%               will only be overwritten within a call to
%               MCMC if Append is set to no.   
%
%     Display - Textual output of algorithm's progress.
%    		[ {on} | off ]
%
%      InDir  - Directory from which initial parameter guesses are obtained. {''}
%    		If empty the initial guess is taken in the function call.
%    		If non-empty the initial guess given in the function call
%    		is ignored and the parameter guesses are taken from the 
%    		file in the directory with the largest label.
%
%    	 MaxT - Maximum number of cycles.  When Blocks is specified 
%    		this is the number of complete cyles through the blocks.
%    		[ pos. integer ]
%
%      Method - Method of generating acceptance probability. 
%    		[ gibbs | m | {mh} ]
%
%      OutDir - Directory in which data are saved. {mcmc_out}
%     
%      Recalc - When 'on' the likelihood, prior etc are recalculated
%               for parameter vectors that were rejected (unchanged)
%               from the previous iteration.  Setting to 'off' can
%               greatly speed up computations.  If, however, the
%               evaluations might change from call to call, for example
%               when data augmentation is being done within a likelihood, then
%               it may be better to set Recalc to 'off'.  It such cases, each
%               iteration typically takes longer to evaluate but the overall 
%               convergence rate is larger due to better mixing properties.
%
%    	Trace - Chart parameter convergence vs. time.
%    		[ {on} | off ]
%
%    WriteOut - Dump Markov chain parameter guesses to a file in OutDir after this
%    		many iterations.
%
%    UseLogLike - The MCMC routine expects that the likelihood function returns the log likelihood.
%                 Use for large data sets.
%      See also MCMC
%

if ((nargout == 0) && (nargin == 0))
	if (nargout == 0)
		disp(['      Append: [ {yes} | no]']);
		disp(['      Blocks: [ pos. integer matrix ]']);
		disp(['DeleteOnStart: [ no | {yes} ]']);
		disp(['     Display: [ {on} | off ]']);
		disp(['       InDir: [ {''} ]']);
		disp(['        MaxT: [ pos. integer ]']);
		disp(['      Method: [ gibbs | m | {mh} ]']);
		disp(['      OutDir: [ {mcmc_out} ]']);
		disp(['      ReCalc: [ {on} | off ]']);
		disp(['       Trace: [ {on} | off ]']);
                disp(['    WriteOut: [ pos. integer ]']);
		disp(['    kdetrack: [ (off) | on ]']);
		disp(['    kdequeue: [ String ]']);
		disp(['    kdesave: [ {mcmc_kde_report} ]']);
        disp(['    UseLogLike: [ {false} | true ]']);
	end
else

	S = struct('Append','yes',...
	           'Blocks',[],...
		   'DeleteOnStart','yes',...
		   'Display','on',...
	           'InDir','',...
		   'MaxT',1000,...
		   'Method','mh',...
		   'OutDir','mcmc_out',...
		   'ReCalc','on',...
		   'Trace','on',...
		   'WriteOut',1000,...
		   'kdetrack','off',...
		   'kdequeue','',...
		   'kdesave','mcmc_kde_report',...
           'UseLogLike',false);
		   
		   
	if (nargin == 0) 
		return;
	else
		if isstruct(varargin{1})
			FieldNames = fieldnames(varargin{1});
			Vals = struct2cell(varargin{1});
		else %A comma-separated list of pairs.
			FieldNames = {varargin{1:2:end-1}};
			Vals = {varargin{2:2:end}};
		end

		for i = 1:length(FieldNames)

			switch lower(FieldNames{i}) 
			case 'append'
				if ~any(strcmp(Vals{i},{'yes','no'}))
					error('Append must be ''yes'' or ''no.''');
				end
			case 'blocks'
				if ~isnumeric(Vals{i})
					error('Blocks must be numeric.');	
				end
				if ~isempty(Vals{i})
					if (~all(all(mod(Vals{i},1) == 0)) | any(any(Vals{i} < 0)))
						error('Block values must non-negative integers.');
					end


					if ~all(diff(sum(Vals{i},1)) == 0)
				        	%This ensures equal coverage...
						%But maybe we don't want this.  Sometimes
						%we might want a subset of the parameters to 
						%be inactive.
						%warning('Block sums not equal for each parameter.')
					end

					for j = 1:size(Vals{i},1)
						ind = find(Vals{i}(j,:)~= 0);
						if ~all(diff(Vals{i}(j,ind)) == 0)
							error('Within rows of Blocks the non-zero entries must be equivalent.');
						end
					end
				end
			case 'deleteonstart'
				if ~any(strcmp(Vals{i},{'yes','no'}))
					error('DeleteOnStart must be ''yes'' or ''no.''');
				end	
			case 'maxt'
				if ~isnumeric(Vals{i})
					error('MaxT must be numeric.');
				end
			case 'method'
				switch lower(Vals{i})
				case 'gibbs'
				case 'm'
				case 'mh'
				otherwise
					error('Unknown method.');
				end
			
			case 'trace'
				if ~ischar(Vals{i})
					error('Trace must be char.');
				end
			case 'display'
				if ~ischar(Vals{i})
					error('Display must be char.');
				end
			case 'indir'
				if ~ischar(Vals{i})
					error('InFile must be character string.')
				end
				
				if ~isempty(Vals{i})
					[DIR,NAME] = fileparts(Vals{i});

					if isempty(DIR)
						DIR = pwd;
					end

					Vals{i} = fullfile(DIR,NAME);
				end
			case 'outdir'
				if ~ischar(Vals{i})
					error('OutFile must be character string.')
				end
				
				if ~isempty(Vals{i})
					[DIR,NAME] = fileparts(Vals{i});
					if isempty(DIR)
						DIR = pwd;
					end
					Vals{i} = fullfile(DIR,NAME);
				end
				
			case 'recalc'
				if ~ischar(Vals{i})
					error('ReCalc must be char.');
				end
			case 'writeout'
				if (~isnumeric(Vals{i}) & (mod(Vals{i},1) == 0))
					error('WriteOut must be integer.')
				end
			case 'kdetrack'
				if ~ischar(Vals{i})
					error('kdetrack must be char.');
				end
			case 'kdequeue'
				if ~ischar(Vals{i})
					error('kdequeue must be char.');
				end
			case 'kdesave'
				if ~ischar(Vals{i})
					error('kdesave must be char.');
                end
            case 'useloglike'
                if ~islogical(Vals{i})
					error('UseLogLike must be a logical.');
				end
			otherwise
				error('Unknown option.');
			end
			F = fieldnames(S);
			
			k = strcmpi(F,FieldNames{i});
			S = setfield(S,F{k},Vals{i});

		end
		
	end	
				
end

