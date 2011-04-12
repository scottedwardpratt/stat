function [] = Configure(FLAG)
%CONFIGURE This function can both set and check for validity that global
%variable CONFIG. CONFIG is used to store all static data necessary to do
%an MCMC search using the related routines. A detailed summary of the
%options are as follows, with each option listing its options, as well as
%default value:
%
%MCMC OPTIONS:
%Append -           Add output file to OutDir during a call to MCMC
%                   If set to 'no' a single output file
%                   will be overwritten at each WriteOut.
%                   Options: true, false    Default: true
%
%AssertSymmetricProposal - If true, the proposal PDF is not calculated. This
%                   options only applies to the Metropolis-Hastings method.
%
%
%Blocks -           Matrix of integers indicating blocking as well as the number 
%                   of iterations taken for each parameter block.
%                   Options:[ pos. integer matrix ]     Default: []
%
%DeleteOnStart -    If 'true', the OutFile (if it exists) is overwritten
%                   at the start of a call to MCMC. The output files  
%                   will only be overwritten within a call to
%                   MCMC if Append is set to 'false'.
%                   Options: true, false    Default: true
%
%Display -          Textual output of algorithm's progress.
%                   Options: true, false    Default: true
%
%InDir  -           Directory from which initial parameter guesses are obtained. {''}
%                   If empty the initial guess is taken in the function call.
%                   If non-empty the initial guess given in the function call
%                   is ignored and the parameter guesses are taken from the 
%                   file in the directory with the largest label.
%                   Options: character string   Default = ''
%
%MaxT -             Maximum number of cycles.  When Blocks is specified 
%                   this is the number of complete cyles through the blocks.
%                   Options: A positive integer     Default = 1000
%
%Method -           Method of generating acceptance probability. 'gibbs' stands for
%                   Gibbs sampling, 'm' for the Metropolis algorithm, 'mh' for the
%                   Metropolis-Hastings algorithm.
%                   Options: 'gibbs','m','mh'   Default: 'mh'
%
%OutDir -           Directory in which data are saved.
%                   Options: character string   Default: 'mcmc_out'
% 
%    
%Recalc -           Guarantees that the likelihood, prior, and  proposal 
%                   likelihood (if Metropolis-Hastings) will be evaluated
%                   for both the candidate and previous parameters.
%                   This is useful under data augmentation, or when stochastic
%                   or other nuissance parameters are desired to be marginalized.
%                   Options: true, false    Default: true
%
%Trace -            Chart parameter convergence vs. time.
%                   Options: true, false    Default: true
%
%WriteOut -         Dump Markov chain parameter guesses to a file in OutDir after this
%                   many iterations.
%                   Options: a positive integer     Default: 100
%
%UseLogLike -       The MCMC routine expects that the likelihood function returns
%                   the log likelihood. Use for large data sets.
%                   Options: true, false        Default: false
%
%GRAPHICS OPTIONS:
%withgraphics-      Determines whether or not an MCMC search shows a
%                   graphical display of its progress.
%                   Options: true, false        Default: true
%
%savegraphics-      Whether or not the displayed graphics are saved.
%                   Options: true, false        Default: true
%
%graphicsfmt-       What format the images are saved in.
%                   Options: 'png','jpg','epsc' Default: 'png'
%
%savegraphicsmod-   Sets a modulus for saving graphics, so that not every
%                   image is saved. For example, setting to 5 would mean
%                   that every 5th image is saved. Increase to decrease
%                   memory used by program.
%
%STATISTICS OPTIONS:
%statsamplesize-    How many iterations taken from the end to run final
%                   statistics over. Warnings given if it is half as large
%                   as total number of iterations, ideally should only
%                   cover the range where the parameter estimates are
%                   relatively constant.
%                   Options: Pos. Int.          Default: 250
%
%DEBUG OPTIONS:
%Note: These options (Debug, Report, Timing) are meant to provide a
%framework for modellers to adapt their own debug/report/timing modes into
%the MADAI frontend. As such, there are few specific built-in capabilities
%to these options. To implement these options, your code should be
%something like:
%
%if(CONFIG.OPTIONS.DEBUG.Debug || CONFIG.OPTIONS.DEBUG.DebugLike)
%DEBUG_MODE = true;
%else
%DEBUG_MODE = false;
%end
%
%
%Debug-             Enters global debug mode. Every debug mode (proposal,
%                   prior, likelihood, factorization, etc.) will be turned
%                   on. Not recommended unless specific debug modes are
%                   unsuccessful.
%                   Options: true, false        Default: false
%
%DebugPrior-        Enables debug mode for prior distribution function.
%                   Options: true, false        Default: false
%
%DebugProp-         Enables debug mode for proposal distribution function.
%                   Options: ture, false        Default: false
%
%DebugLikelihood-   Enables debug mode for likelihood distribtuion
%                   function.
%                   Options: true, false        Default: false
%
%DebugPost-         Enables debug mode for posterior distrubtion function.
%                   Options: true, false        Default: false
%
%DebugFactor-       Enables debug mode for cholesky factorization function
%                   for cosmosurvey.
%                   Options: true, false        Default: false
%
%ReportLike-        Enables report mode for likelihood distribution
%                   function.
%                   Options: true, false        Default: false
%ReportPrior-       Enables report mode for prior distribution function.
%                   Options: true, false        Default: false
%
%ReportFactor-      Enables report mode for cholesky factorization for
%                   cosmosurvey project.
%                   Options: true, false        Default: false
%
%ReportPost-        Enables report mode for posterior distribution
%                   function.
%                   Options: true, false        Default: false
%
%ReportProp-        Enables report mode for proposal factorization
%                   function.
%                   Options: true, false        Default: false
%
%ReportPrior-       Enables report mode for prior distribution function.
%                   Options: true, false        Default: false
%
%HandleMessages-    Enables support of the HandleMessages m-file. This
%                   allows debug messages and report messages to be
%                   displayed to the screen and stored in the global 
%                   variable LOGS so that they can be uploaded with
%                   the trace. See HandleMessages.m for more details.
%                   Options: true,false         Default: true
%
%LogDebug-          Determines if HandleMessages logs debug messages.
%                   Options: true,false         Default: true
%
%LogMessages-       Determiens if all other messages are logged by
%                   HandleMessages.
%                   Options: true,false         Default: true
%
%TIMING OPTIONS:
%Timing-            Enables global timing mode for all functions. Not
%                   recommended unless all timing data is needed.
%                   Options: true, false        Default: false
%
%TimeLike-          Enables timing mode for likelihood distribution.
%                   Options: true, false        Default: false
%
%TimeFactor-        Enables timing mode for cholesky factorization.
%                   Options: true, false        Default: false
%
%TimeProp-          Enables timing mode for proposal distribution.
%                   Options: true, false        Default: false
%
%TimePrior-         Enables timing mode for prior distribution.
%                   Options: true, false        Default: false
%
%TimeIteration-     Enables a mode for reporting the time it takes to go
%                   through one complete iteration.
%                   Options: true, false        Default: false
%
%LIKELIHOOD OPTIONS:
%DetMethod-         Determines which method the likelihood function uses to
%                   calculate the determinant of the covariance matrix.
%                   Leaving this field blank uses the brute force, sum of
%                   the diagonals approach. 'Sampled' uses the sampling
%                   approach from Dr. Dougherty, adapting the Pace & LaSage
%                   approach. 'novak' is a similar method with different code
%                   optimizing it for small matrices.
%                   Options: 'Sampled' 'novak' ''   Default: ''
%
%Qsize-             Determines how big the queue of stored recent
%                   factorizations is. Recommended to be at least 2,
%                   increasing this causes more memory to be used.
%                   Options: pos. int.              Default: 3
%
%Cusping-           Whether of not the covariance matrix accomodates
%                   cusping caused by spatial data. In the cosmosurvey
%                   project, the cusp locations are determined by the nl
%                   value.
%                   Options: true, false            Default: false
%
%SampleBlockSize-   How large the sample blocks are in the sampled
%                   determinant.
%                   Options: true, false            Default: 20
%
%SampleOffset-      How many initial rows of the covariance matrix are 
%                   skipped when calculating the sampled determinant. Set
%                   this to 1 to not skip any rows.
%                   Options: pos. int.              Default: 1
%See also MCMC

global CONFIG;
global LOGS;

switch FLAG
    case 'set' %set the parameters to their default values. See above to determine exactly what each option does
        CONFIG.OPTIONS = [];
        CONFIG.LOGPARAM = [];
        CONFIG.STATIC = [];
        CONFIG.CHOL = [];
        CONFIG.META = [];
        
        CONFIG.META = struct('mcmcrunnickname','cosmo_run',...
                             'runcomment','',...
                             'snifflibpath','\Users\kevinnovak\',...
                             'snifflib_version',1.1,...
                             'trace_config','',...
                             'param_config','',...
                             'trace_repository','',...
                             'param_repository','');
                            
        CONFIG.OPTIONS.MCMC = struct('Append',true,...
                                     'Blocks',[],...
                                     'DeleteOnStart',true,...
                                     'Display',true,...
                                     'InDir','',...
                                     'MaxT',1000,...
                                     'Method','mh',...
                                     'OutDir','mcmc_out',...
                                     'ReCalc',true,...
                                     'Trace',true,...
                                     'WriteOut',100,...
                                     'kdetrack',false,...
                                     'kdequeue','',...
                                     'AssertSymmetricProposal',false,...
                                     'kdesave','mcmc_kde_report',...
                                     'UseLogLike',false);
        
        CONFIG.OPTIONS.GRAPHICS = struct('withgraphics',true,...
                                         'savegraphics',true,...
                                         'graphicsfmt','png',...
                                         'savegraphicsmod',1,...
                                         'DeleteGraphicsOnStart',true);
        CONFIG.OPTIONS.STATS = struct('statsamplesize',250);
        CONFIG.OPTIONS.DEBUG = struct('Debug',false,...
                                      'DebugPrior',false,...
                                      'DebugProp',false,...
                                      'DebugLikelihood',false,...
                                      'DebugPost',false,...
                                      'DebugFactor',false,...
                                      'ReportLike',false,...
                                      'ReportFactor',false,...
                                      'ReportPost',false,...
                                      'ReportProp',false,...
                                      'ReportPrior',false,...
                                      'HandleMessages',true,...
                                      'LogDebug',true,...
                                      'LogMessages',true,...
                                      'LogErrors',true,...
                                      'LogWarnings',true);
        CONFIG.OPTIONS.TIMING = struct('Timing',false, ...
                                       'TimeLike',false,...
                                       'TimeFactor',false,...
                                       'TimeProp',false,...
                                       'TimePrior',false,...
                                       'Timeiteration',false);
        CONFIG.OPTIONS.LIKELIHOOD = struct('DetMethod','',...
                                           'Qsize',3,...
                                           'Cusping',false,...
                                           'SampleBlockSize',20,...
                                           'SampleOffset',1);
        CONFIG.OPTIONS.COSMO = struct('savetrace',true);
        CONFIG.PARAMETERS = {};
        CONFIG.FUNCTIONS = struct('Likelihood','',...
                                  'Proposal','',...
                                  'Prior','',...
                                  'Posterior','',...
                                  'PreMCMC','',...
                                  'PostMCMC','');
                              
        CONFIG.LIKELIHOOD = struct('DATAFILE', '');
        CONFIG.PROPOSAL = struct('MIXING_STDDEV', '');
        

        
    case 'check' %check to make sure options are valid
        
        LOGS.Debug = cell(CONFIG.OPTIONS.MCMC.MaxT, 1);
        LOGS.Messages = cell(CONFIG.OPTIONS.MCMC.MaxT,1);
        LOGS.Error = '';
        
        %Known bug: If this option is entered incorrectly, the generated
        %error message will not be stored in the LOGS variable.
        if(CONFIG.OPTIONS.DEBUG.HandleMessages)
            HM = true;
        else
            HM = false;
        end
        
        CONFIG_NAMES = fieldnames(CONFIG);
        
        for a = 1:length(CONFIG_NAMES)
            
            switch(CONFIG_NAMES{a})
                case 'OPTIONS'
                    F = fieldnames(CONFIG.OPTIONS);

                    for j = 1:length(F)
                        if(any(strcmp((F{j}), { 'MCMC' , 'GRAPHICS' , 'STATS' , 'DEBUG' , 'TIMING' , 'COSMO','LIKELIHOOD'})) == 1)
                            FieldNames = fieldnames(CONFIG.OPTIONS.(F{j}));
                            Vals = struct2cell(CONFIG.OPTIONS.(F{j}));

                            for i = 1:length(FieldNames)

                                switch lower(FieldNames{i}) 
                                    %True/False options
                                    case {'append', 'debug', 'debugprior', 'debugpost','debuglikelihood', 'debugprop',...
                                          'reportlike', 'reportprop', 'reportprior', 'reportpost','deleteonstart', ...
                                          'trace', 'display','recalc','kdetrack', 'useloglike', 'withgraphics',...
                                          'timing', 'timelike', 'timefactor', 'timeprop','timeprior',...
                                          'timeiteration', 'savegraphics','cusping','debugfactor','reportfactor',...
                                          'handlemessages','savetrace','logdebug','logmessages','deletegraphicsonstart',...
                                          'assertsymmetricproposal', 'logerrors','logwarnings'}
                                        if ~islogical(Vals{i})
                                            errormsg = [FieldNames{i} ' must be a logical.'];
                                            if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                            error(errormsg);
                                        end
                                    %positive integer options
                                    case {'maxt', 'writeout','statsamplesize', 'savegraphicsmod','qsize',...
                                          'sampleblocksize','sampleoffset'}
                                        if(mod(Vals{i},1) ~= 0 || Vals{i} < 1)
                                            errormsg = [FieldNames{i} ' must be a positive integer, it is: ' num2str(Vals{i}) '.'];
                                            if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                            error(errormsg);
                                        end
                                    %special for blocks options    
                                    case 'blocks'
                                        if (~isempty(Vals{i})&&~isinteger(Vals{i}))
                                            errormsg = [FieldNames{i} ' must be numeric.'];
                                            if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                            error(errormsg);
                                        end
                                        if ~isempty(Vals{i})
                                            if (~all(all(mod(Vals{i},1) == 0)) || any(any(Vals{i} < 0)))
                                                errormsg = 'Block values must non-negative integers.';
                                                if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                error(errormsg);
                                            end


                                            if ~all(diff(sum(Vals{i},1)) == 0)
                                                %This ensures equal coverage...
                                                %But maybe we don't want this.  Sometimes
                                                %we might want a subset of the parameters to 
                                                %be inactive.
                                                %warning('Block sums not equal for each parameter.')
                                            end

                                            for k = 1:size(Vals{i},1)
                                                ind = Vals{i}(k,:)~= 0;
                                                if ~all(diff(Vals{i}(k,ind)) == 0)
                                                    errormsg = 'Within rows of Blocks the non-zero entries must be equivalent.';
                                                    if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                    error(errormsg);
                                                end
                                            end
                                        end
                                    %options with a few valid choices
                                    case 'method'
                                        switch lower(Vals{i})
                                            case 'gibbs'
                                            case 'm'
                                            case 'mh'
                                            otherwise
                                                errormsg = 'Unknown mcmc method.';
                                                if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                error(errormsg);
                                        end
                                    case 'graphicsfmt'
                                        switch lower(Vals{i})
                                            case 'png'
                                            case 'jpg'
                                            case 'epsc'
                                            otherwise
                                                errormsg = 'Unknown graphics format.';
                                                if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                error(errormsg);
                                        end
                                    case 'detmethod'
                                        switch lower(Vals{i})
                                            case 'sampled'
                                            case 'novak'
                                            case ''
                                            otherwise
                                                errormsg = 'Unknown determinant method.';
                                                if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                error(errormsg);
                                        end
                                    case 'propmethod'
                                        switch lower(Vals{i})
                                            case 'rand'
                                            case 'pdf'
                                            otherwise
                                                errormsg = 'Unknown proposal method.';
                                                if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                error(errormsg);
                                        end
                                    %file directories
                                    case {'indir', 'outdir'}
                                        if ~ischar(Vals{i})
                                            errormsg = [FieldNames{i} ' must be character string.'];
                                            if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                            error(errormsg);
                                        end

                                        if ~isempty(Vals{i})
                                            [DIR,NAME] = fileparts(Vals{i});

                                            if isempty(DIR)
                                                DIR = pwd;
                                            end

                                            Vals{i} = fullfile(DIR,NAME);
                                        end
                                    %strings
                                    case {'kdequeue', 'kdesave'}
                                        if ~ischar(Vals{i})
                                            errormsg = [FieldNames{i} ' must be char.'];
                                            if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                            error(errormsg);
                                        end
                                    %structures
                                    case 'mixing_stddev'
                                        if ~isempty(Vals{i})
                                            if ~isstruct(Vals{i} && ~iscell(Vals{i}))
                                                errormsg = [FieldNames{i} ' must be structure.'];
                                                if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                                error(errormsg);
                                            end
                                        end
                                    case 'parameters'
                                        if(~iscell(Vals{i}) || isempty(Vals{i}))
                                            errormsg = [FieldNames{i} ' must be a cell array of values'];
                                            if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                            error(errormsg);
                                        end
                                    otherwise
                                        errormsg = ['Unknown option ' FieldNames{i} '.'];
                                        if(HM), HandleMessages(3, errormsg, ['-' FieldNames{i} ]); end;
                                        error(errormsg);
                                end
                            end
                        else
                            errormsg = ['Unknown options set ' F{j}];
                            if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                            error(errormsg);
                        end
                    end
                case {'CHOL','LOGPARAM'}
                    if(~isempty(CONFIG.(CONFIG_NAMES{a}))&&~iscell(CONFIG.(CONFIG_NAMES{a})))
                        errormsg = [CONFIG_NAMES{a} ' needs to be a cell array of values.'];
                        if(HM), HandleMessages(3, errormsg, ['-' CONFIG_NAMES{a} ]); end;
                        error(errormsg);
                    end
                case {'STATIC','LIKELIHOOD','PROPOSAL'}
                    if(~isstruct(CONFIG.(CONFIG_NAMES{a})))
                        errormsg = [CONFIG_NAMES{a} ' needs to be a structure of options and values.'];
                        if(HM), HandleMessages(3, errormsg, ['-' CONFIG_NAMES{a} ]); end;
                        error(errormsg);
                    end
                case 'FUNCTIONS'
                    F = fieldnames(CONFIG.FUNCTIONS);
                    Vals = struct2cell(CONFIG.FUNCTIONS);
                    
                    for i = 1:length(F)
                        switch F{i}
                            case{'Likelihood','Proposal','Prior'}
                                if(~ischar(Vals{i})||isempty(Vals{i}))
                                    errormsg = [F{i} ' must not be empty and must be a character string.'];
                                    if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                    error(errormsg);
                                end
                            case{ 'PreMCMC','PostMCMC','Posterior'}
                                if(~isempty(Vals{i})&&~ischar(Vals{i}))
                                    errormsg = [F{i} ' must be a character string.'];
                                    if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                    error(errormsg);
                                end
                            otherwise
                                errormsg = ['Unknown field ' F{i}];
                                if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                error(errormsg);
                        end
                    end
                    
                case 'META'
                    F = fieldnames(CONFIG.META);
                    Vals = struct2cell(CONFIG.META);
                    
                    for i = 1:length(F)
                        switch F{i}
                            case {'mcmcrunnickname', 'snifflibpath'}
                                if(~ischar(Vals{i}) || isempty(Vals{i}))
                                    errormsg = [F{i} ' must not be empty and must be a character string.'];
                                    if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                    error(errormsg);
                                end
                            case 'runcomment'
                                if(~ischar(Vals{i}))
                                    errormsg = [F{i} ' must be a character string'];
                                    if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                    error(errormsg);
                                end
                            case 'snifflib_version'
                                if(~isfloat(Vals{i}))
                                    errormsg = [F{i} ' must be a floating point number.'];
                                    if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                    error(errormsg);
                                end
                            case {'trace_config','param_config'}
                                switch Vals{i}
                                    case 'xml'
                                    case 'sql'
                                    case ''
                                    otherwise
                                        errormsg = ['Unknown configuration in ' F{i}];
                                        if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                        error(errormsg);
                                end
                            case {'trace_repository','param_repository'}
                                if(~isempty(Vals{i}) && ~ischar(Vals{i}))
                                    errormsg = [ F{i} ' must either be empty or a character string.'];
                                    if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                    error(errormsg);
                                end
                            otherwise
                                errormsg = ['unknown option ' F{i} ];
                                if(HM), HandleMessages(3, errormsg, ['-' F{i} ]); end;
                                error(errormsg);
                        end
                    end
                case 'PARAMETERS'
                    if(~iscell(CONFIG.PARAMETERS))
                        errormsg = 'CONFIG.PARAMETERS must be a cell array of strings.';
                        if(HM), HandleMessages(3, errormsg, '-PARAMETERS'); end;
                        error(errormsg);
                    end
                otherwise
                    errormsg = ['Unknown CONFIG field ' CONFIG_NAMES{a}];
                    if(HM), HandleMessages(3, errormsg, ['-' CONFIG_NAMES{a} ]); end;
                    error(errormsg);
            end
        end
    otherwise
        errormsg = ['Unknown flag ' FLAG];
        if(CONFIG.OPTIONS.DEBUG.HandleMessages), HandleMessages(3, errormsg, ['-' FLAG ]); end;
        error(errormsg);
end
end
