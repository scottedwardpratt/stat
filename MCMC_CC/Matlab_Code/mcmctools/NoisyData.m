function [ ] = NoisyData(filename, nl, nz, alpha, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global CONFIG;

fid = fopen(filename);
if(fid == -1)
    if(CONFIG.OPTIONS.DEBUG.HandleMessages)
        HandleMessages(3, 'NoisyData: Cannot open data file.');
    end
    error('NoisyData: Cannot open data file.');
end

%Set MU equal to model output. Most likely, additional work will be needed
%here to convert from a string to a vector.
MU = cell2mat(textscan(fid,'%.15f'));
frewind(fid);
%Form a spatial grid to add noise (assuming spatially dependant noise)
xx = linspace(-3,3,nz);
yy = linspace(-3,3,nl);
[X,Y] = ndgrid(xx,yy);

D = squareform(pdist([X(:),Y(:)],'euclidean'));

%Form factorization for noise based on MESS method
Salphainv = expm(-alpha*D);
VARCOV = Salphainv'*Salphainv;
VARCOV = ((sigma^2)./VARCOV(1,1)).*VARCOV;

%Use MATLAB's mvnrnd in stats toolbox to handle issues of near-negative
%eigen values due to numerical errors...
Zinv = mvnrnd(MU(:),VARCOV);

fwrite(fid, Zinv);
fclose(fid);
end

