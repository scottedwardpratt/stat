function [ OUT, AUG ] = cosmoLikeLihood(THETA,AUG)
%LikeLihood Evaluates the Likelihood function at a given point
%   Depends on a function called BuildCholeskyFinal to recusively build a
%   physical toeplitz covariance matrix, then uses a MVN likelihood
%   function to evaluate the log likelihood at a given point in parameter
%   space.

%************************PARAMETERS************************
%MCMC parameters:
persistent CALLCOUNT;
global CONFIG;

%Determine Cosmoplot Parameters:
COSMOPARAMS = setdiff(fieldnames(THETA),CONFIG.CHOL);

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

%************************END PARAMETERS ********************
if(DEBUG)
    if(HM)
        HandleMessages(1,'Likelihood called.');
    else
        disp('Likelihood called.');
    end
end

N= CONFIG.STATIC.nl*CONFIG.STATIC.nz; %Total number of data points

if (isempty(CALLCOUNT))
	CALLCOUNT = 0;
end
CALLCOUNT = CALLCOUNT + 1;

for i=1:length(CONFIG.LOGPARAM)
	THETA.(CONFIG.LOGPARAM{i}) = exp(THETA.(CONFIG.LOGPARAM{i}));
end

%************************CHOLESKY FACTORIZATION*************
%Create Cholesky factorization first, saving time by making sure parameters
%are valid before running the model.
if(TIMING)
    tic;
end

if(DEBUG)
    if(HM)
        HandleMessages(1, 'Getting Factorization');
    else
        disp('Getting Factorization');
    end
end

V = zeros(1,N);

for k = 1:length(CONFIG.CHOL)
    if(mod(str2num(CONFIG.CHOL{k}(2)),1)~= 0)
        if(HM), HandleMessages(3, 'Error interpreting cholesky parameters.', '-x*1'); end;
        error('Error interpreting cholesky parameters.');
    end
    index = str2num(CONFIG.CHOL{k}(2));

    V(index) = THETA.(CONFIG.CHOL{k});
end

L = sparse(BuildCholeskyFinal(V));

if( (any(isnan(L)) == 1))
    if(HM), HandleMessages(3, 'NaNs detected in factorization.'); end;
    error('NaNs detected in factorization');
end
if(any(isinf(L))==1)
    if(HM), HandleMessages(3, 'Infs detected in factorization'); end;
    error('Infs detected in factorization');
end
if(DEBUG)
    if(HM)
        HandleMessages(1, 'Done.');
    else
        disp('Done.');
    end
end

if(TIMING)
    times(1)=toc;
end
%***********************RUNNING MODEL***********************
if(DEBUG)
    if(HM)
        HandleMessages(1, 'Querying Model...');
    else
        disp('Querying Model...');
    end
end
if(TIMING)
    tic;
end

F_STATIC = fieldnames(CONFIG.STATIC);

CMD = '';
%Add variable parameters to system command
for i=1:length(COSMOPARAMS)
	CMD = [CMD,' -',COSMOPARAMS{i},' ', num2str(THETA.(COSMOPARAMS{i}), '%.15f')];
end
%Add static parameters to system command
for j = 1:length(F_STATIC)
    if(any(strcmp(F_STATIC{j},{'nz','nl'})) == 1)
        CMD = [CMD, ' -' F_STATIC{j} ' ' num2str(CONFIG.STATIC.(F_STATIC{j}), '%f')];
    else
        CMD = [CMD, ' -' F_STATIC{j} ' ' num2str(CONFIG.STATIC.(F_STATIC{j}), '%.15f')];
    end
end
 
CMD = ['/usr/local/bin/cosmosurvey ',CMD];

[status,result] = system(CMD);

if(status ~=0)
    if(HM), HandleMessages(3, 'Error in running model.'); end;
    error('Error in running model.');
end


if(DEBUG)
    if(HM)
        HandleMessages(1, 'Done.');
    else
        disp('Done.');
    end
end

modeldata = cell2mat(textscan(result,'%f'));

if(any(modeldata==0) == 1)
    modeldata = log(modeldata + 0.5);
else
    modeldata = log(modeldata);
end

% if(any(isnan(modeldata))==1 || any(isinf(modeldata))==1)
%     error('NaNs  or Infs detected in model data.');
% end

if(DEBUG)
    if(HM)
        HandleMessages(1, 'Data converted to numbers.');
    else
        disp('Data converted to numbers.');
    end
end
if(TIMING)
    times(2)=toc;
end
%*********************GET LOG DETERMINANT*********************
if(TIMING)
    tic;
end

switch lower(CONFIG.OPTIONS.LIKELIHOOD.DetMethod)
    
    case 'sampled'
        if(DEBUG)
            if(HM)
                HandleMessages(1, 'Using sampled det.');
            else
                disp('Using sampled det.');
            end
        end
        SampleLogDet = [];
        start = CONFIG.OPTIONS.LIKELIHOOD.SampleOffset;
        k=1;
        while (start<=((N-CONFIG.OPTIONS.LIKELIHOOD.SampleBlockSize)+1))
        	rows = start:start+CONFIG.OPTIONS.LIKELIHOOD.SampleBlockSize-1;
        	sample = sub2ind(size(L),rows,rows);
        	start = start + 2*CONFIG.OPTIONS.LIKELIHOOD.SampleBlockSize;
        	%We jump 2 blocks so as to increase/ensure independence of our samples...
        	SampleLogDet(k) = sum(log(abs(L(sample))))/CONFIG.OPTIONS.LIKELIHOOD.SampleBlockSize;
            if(imag(SampleLogDet(k)) > 0)
                if(HM), HandleMessages(3, 'imaginary determinant'); end;
                error('imaginary determinant');
            end
        	k = k+1;
        end
        
        %Taking the "average of the averages" is a hack to avoid doing lots of tests about
        %independence and variance checking and so on.
        
        BlockAvgLogDet = mean(SampleLogDet);
        
        LogDet = 2*N*BlockAvgLogDet;
    otherwise
        if(DEBUG) 
            if(HM)
                HandleMessages(1, 'Using Brute Force Method.');
            else
                disp('Using Brute Force Method');
            end
        end
        LogDet = 2*sum(log(abs(diag(L))));
end

if(DEBUG)
    if(HM)
        HandleMessages(1, 'Determinant Calculated.');
    else
        disp('Determinant Calculated');
    end
end

if(TIMING)
    times(3) = toc;
end

%*******************CALCULATE EXPONENTIAL *********************
if(TIMING)
    tic;
end

%open data file
fid = fopen(CONFIG.LIKELIHOOD.DATAFILE);
if(fid == -1)
    if(HM), HandleMessages(3, 'Cannot open data file.'); end;
    error('Cannot open data file.');
end
%read in numbers, do some massaging to get them readable

data= cell2mat(textscan(fid,'%f'));

if(any(data == 0) == 1)
    data = log(data + 0.5);
else
    data = log(data);
end

if(any(isnan(data))==1)
    if(HM), HandleMessages(3, 'NaNs detected in synthetic data.'); end;
    error('NaNs detected in synthetic data.');
end

fclose(fid);

if(DEBUG)
    if(HM)
        HandleMessages(1, 'Data read in.');
    else
        disp('Data read in.');
    end
end

if(size(modeldata,1)~=size(data,1))
    if(HM), HandleMessages(3, ['Prediction has size ' num2str(size(modeldata)) ' and data has size ' num2str(size(data))]); end;
    error(['Prediction has size ' num2str(size(modeldata)) ' and data has size ' num2str(size(data))]);
end
LinvX = L\(data-modeldata);

Exponent = LinvX'*LinvX;

if(any(isnan(LinvX))==1)
    if(HM), HandleMessages(3, 'NaNs detected in exponent'); end;
    error('NaNs detected in exponent');
end

if(DEBUG)
    if(HM)
        HandleMessages(1, 'Exponent Calculated.');
    else
        disp('Exponent Calculated.');
    end;
end
if(TIMING)
    times(4) = toc;
end
%******************PUTTING IT ALL TOGETHER ********************
if(TIMING)
    tic;
end

zz = linspace(CONFIG.STATIC.z1, CONFIG.STATIC.z2,CONFIG.STATIC.nz);
ll = linspace(CONFIG.STATIC.l1, CONFIG.STATIC.l2,CONFIG.STATIC.nl);
[Z,LUM] = meshgrid(zz,ll);


if (CONFIG.OPTIONS.GRAPHICS.withgraphics)
    if(mod(CALLCOUNT, CONFIG.OPTIONS.GRAPHICS.savegraphicsmod)==0||CALLCOUNT == 1)
        dataplot = reshape(data,CONFIG.STATIC.nl,CONFIG.STATIC.nz);
        modelplot = reshape(modeldata,CONFIG.STATIC.nl, CONFIG.STATIC.nz);
        figure(2);
       
        %zsize = size(Z)
        %lsize = size(LUM)
        %datasize = size(dataplot)
        
        hold off;
        H1 = surf(Z,LUM,dataplot);
        hold on;
        H2 = surf(Z,LUM,modelplot);
        alpha(H1,0.5);
        alpha(H2,0.5);
        xlabel('Redshift');
        ylabel('Luminosity');
        zlabel('Log(Counts)');
        axis tight;
        hold off;
        colormap gray;
        colorbar;
        drawnow;
        if (CONFIG.OPTIONS.GRAPHICS.savegraphics)
            FILE = fullfile(CONFIG.OPTIONS.MCMC.OutDir,[CONFIG.META.mcmcrunnickname '_',sprintf('%d',CALLCOUNT)]);
            if(DEBUG)
                if(HM)
                    HandleMessages(3, ['Storing ' CONFIG.META.mcmcrunnickname '_' sprintf('%d',CALLCOUNT)]);
                else
                    disp(['Storing ' 'cosmo_' sprintf('%d',CALLCOUNT)]);
                end
            end
            switch CONFIG.OPTIONS.GRAPHICS.graphicsfmt
                case 'png'
                    print('-dpng',FILE);
                case 'jpg'
                    print('-djpeg',FILE);
                case 'epsc'
                    print('-depsc',FILE); 
                otherwise
                    if(HM)
                        HandleMessages(4, ['Graphics format ' CONFIG.OPTIONS.GRAPHICS.graphicsfmt ' not recognize. Not saved.']);
                    end
                    warning(['Graphics format ' CONFIG.OPTIONS.GRAPHICS.graphicsfmt ' not recognize. Not saved.']);
            end
        end
    end
else
    CALLCOUNT = CALLCOUNT-1;
end

OUT = -(N/2)*log(2*pi) - (0.5)*LogDet -(0.5)*Exponent;

if(~CONFIG.OPTIONS.MCMC.UseLogLike)
    OUT = exp(OUT);
end

if(TIMING)
    times(5)=toc;
end

if(CONFIG.OPTIONS.DEBUG.ReportLike)
    if(HM)
        HandleMessages(2, ['Reporting Likelihood Parts:\nLogDet: ' num2str(LogDet) '\nExponent: ' num2str(Exponent) '\nOUT: ' num2str(OUT)]);
    else
        sprintf(['Reporting Likelihood Parts:\nLogDet: ' num2str(LogDet) '\nExponent: ' num2str(Exponent) '\nOUT: ' num2str(OUT)]);
    end
end

if(DEBUG)
    if(HM)
        HandleMessages(1, ['Finished, likelihood is ' num2str(OUT)]);
    else
        disp(['Finished, likelihood is ' num2str(OUT)]);
    end
end
if(TIMING)
    if(HM), HandleMessages(2, ['Factorization took: ' num2str(times(1)) ' seconds\n' 'Running model took: ' num2str(times(2)) ' seconds\n' 'Calculating determinant took: ' num2str(times(3)) ' seconds\n' 'Calculating exponential took: ' num2str(times(4)) ' seconds\n' 'Saving graphic and putting it together: ' num2str(times(5)) ' seconds']);end;
    disp(['Factorization took: ' num2str(times(1)) ' seconds']);
    disp(['Running model took: ' num2str(times(2)) ' seconds']);
    disp(['Calculating determinant took: ' num2str(times(3)) ' seconds']);
    disp(['Calculating exponential took: ' num2str(times(4)) ' seconds']);
    disp(['Saving graphic and putting it together: ' num2str(times(5)) ' seconds']);
end

