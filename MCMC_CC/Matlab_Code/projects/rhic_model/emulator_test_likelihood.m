function [like, AUG] = emulator_test_likelihood(THETA, AUG)

global CONFIG;
global DATA;

if(CONFIG.OPTIONS.TIMING.Timing || CONFIG.OPTIONS.TIMING.TimeLike)
    TIMING = true;
else
    TIMING = false;
end

if(CONFIG.OPTIONS.DEBUG.Debug || CONFIG.OPTIONS.DEBUG.DebugLikelihood)
    DEBUG = true;
else
    DEBUG = false;
end
if(CONFIG.OPTIONS.DEBUG.HandleMessages)
    HM = true;
else
    HM = false;
end

if(TIMING)
    times = zeros(1,5);
end

%***********************RUNNING MODEL***********************
if(DEBUG)
    if(HM)
        HandleMessages(1, 'Querying Emulator...');
    else
        disp('Querying Emulator...');
    end
end
if(TIMING)
    tic;
end

[meanvals, errors] = QueryEmulator([getenv('SCRIPT_HOME') '/../model-data/rhic-2nd-test'], THETA);

if(DEBUG)
    if(HM)
        HandleMessages(1, 'Emulator output converted to numbers.');
    else
        disp('Emulator output converted to numbers.');
    end
end
if(TIMING)
    times(1)=toc;
end

%*************READING IN DATA********************
if(DEBUG)
    if(HM)
        HandleMessages(1, 'Reading In Data...');
    else
        disp('Reading In Data...');
    end
end
if(TIMING)
    tic;
end

data = DATA;

SIGMA = diag(abs(errors));

if(TIMING)
    times(2)=toc;
end
%*************PUTTING IT TOGETHER****************
if(DEBUG)
    if(HM)
        HandleMessages(1, 'Calculating Final Likelihood...');
    else
        disp('Calculating Final Likelihood...');
    end
end
if(TIMING)
    tic;
end

like = mvnpdf(data-meanvals, 0, SIGMA);


if(CONFIG.OPTIONS.MCMC.UseLogLike)
    like = log(like);
end

if(TIMING)
    times(3)=toc;
end

if(DEBUG)
    if(HM)
        HandleMessages(1, ['Finished, likelihood is ' num2str(like)]);
    else
        disp(['Finished, likelihood is ' num2str(like)]);
    end
end
if(TIMING)
    if(HM), HandleMessages(2, ['Emulator took: ' num2str(times(1)) ' seconds\n' 'Reading in data took: ' num2str(times(2)) ' seconds\n' 'Calculating distribution took: ' num2str(times(3)) ' seconds']);end;
    disp(['Emulator took: ' num2str(times(1)) ' seconds']);
    disp(['Reading in data took: ' num2str(times(2)) ' seconds']);
    disp(['Calculating distribution took: ' num2str(times(3)) ' seconds']);
end

end

%add in hierarchical errors.
