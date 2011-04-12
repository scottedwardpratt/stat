function OUT = warehousefun(DIR,C,I,FUN,varargin)
%WAREHOUSEFUN Calculate a function over data stored in a warehouse directory.
%   WAREHOUSEFUN(DIRECTORY,C,I,FUN,P1,P2,...)
%  
%   DIRECTORY the path of the data warehouse.  This directory should
%   contain matfiles with names following the convention 
%   out_i_j_k.mat  where i indicates the sequence with which the file 
%   was created, j is the index of the first observation in the 
%   file and k is the index of the last observation in the file.  Thus
%   'out_7_97_302' indicates that the seventh output which contains observations
%   97 through 302. 
%
%   Subsetting parameters C and I can be used to obtain various effects.  
%   Both parameters may be left empty to achieve the default effect of including
%   all observations.  Other effects are the following: 
%      Cut-off limit C:
%         positive integer - observations >= C are included
%         negative integer - observations <= C are included.
%                    empty - no cutoff limits are used
%      Subset vector I:
%                    empty - all observations considered
%           integer vector - only consider these observations 
%   
%   FUN should having calling syntax FUN(X,t,FLAG,P1,P2,...)
%   where X is an array of data with one row per observation
%   and t holds the observation indices.  
%   Note that X may be a matrix, an array of structs
%   or a cell array. 
%
%   FLAG is 'start' when the
%   first file is read, '' while intermediate files are being read 
%   FLAG = 'report' when the last file is read.  The value of
%   FLAG can be used to determine when to carry out certain calculations
%   or may be ignored.
%   
%  


%This should be stream lined for the case when all data in DIR is to be read.
%That will be the most common case and we don't want to have to keep doing
%unnceccary ckecking.
%

S = dir(fullfile(DIR,'out_*.mat'));

t = 0;
N = 0;

FILTER = '';

if (nargin >= 4)
	FILTER = varargin{1};
	PARAMS = {varargin{2:end}};
end

FLAG = 'start';

if (length(S) == 0)
	error(['Can''t get guesses -- no OUT files in',DIR,'.']);
end


%Note that the convention for the mat files in the data warehouse
%are as follows
%out_<positive nonzero integer sequence number>_<minimum observation>_<maximum observation>
%
%For example: 
%   out_3_51_103.mat  implies that this mat file is the 3rd set of ouput (data)
%   and contains ouputs (observations) 51 through 103.   
%
%

NUM = parsenum({S.name});%Get the number information from file names.
NUM = [NUM,[1:size(NUM,1)]']; %Stick an index on in the last column which will be into the list of file
				%names

if (length(unique(NUM(:,1)))< size(NUM,1))
	error(['The sequence numbers of the ouput files in ''',DIR,''' are not unique.']);
end

ALLOBS = logical(0);

if isempty(I) %Take I to be all observations.
	ALLOBS = logical(1); %True regardless if there is a cutoff.
	I = [min(NUM(:,2)):max(NUM(:,3))];
elseif isequal(I,[min(NUM(:,2)):max(NUM(:,3))])
	%The use input all observations manually. 
	ALLOBS = logical(1); %True regardless if there is a cutoff.
end


if ~isempty(C) %There is a cut-off.  We use this to limit the files that we are considering outright.
	
	if ~isequal(size(C),[1 1])
		error('Cutoff must be a scalar.');
	end
	
	if (C > 0)
		if ALLOBS
			%Keep files where C is less than every observation or 
			%where C is in the middle somewhere.
			NUM = NUM((NUM(:,3) >=C),:);
		end
		I = I(I >= C);
	else
		
		
		if ALLOBS
			%Keep files where C is greater than every observation or 
			%where C is in the middle somewhere.
			NUM = NUM((NUM(:,2) <=abs(C)),:);
		end
		I = I(I <= abs(C));
	end
end


%Now we consider observations within the remaining files...
for i = 1:size(NUM,1) %Iterate over the out_*.mat files.

	
	if ALLOBS
		OPEN = logical(1);
	elseif any((I>=NUM(i,2)) & (I<=NUM(i,3)))
		OPEN = logical(1);
	else
		OPEN = logical(0);
	end
		
	if OPEN %Open this file. Now we need to pick out the observations that we should.	
		
		if (i == size(NUM,1))
			FLAG = 'report';
		elseif (i == 1)
			FLAG = 'start';
		else
			FLAG = '';
		end
		

		X = getfield(load(fullfile(DIR,S(NUM(i,4)).name),'X'),'X');
		t = [NUM(i,2):NUM(i,3)];
		
		
		ind = logical(ones(1,size(X,1))); %The default is to use all observations.
		
		if ~isempty(C) %Deal with cutoff.
			if (C<0)
				ind = (t <= C);
			else
				ind = (t >= C);
			end
			ind(ind) = ismember(t(ind),I);
		
		else %C is empty.
			if ~ALLOBS
				ind(ind) = ismember(t(ind),I);
			end	
		end
			
		
		
		OUT = feval(FUN,X(ind,:),t(ind),FILTER,FLAG,PARAMS{:});
	end
end

clear(FUN); %It will often be the case that the user will put persistent variables in the FUN
	%to calculate running average or running standard deviations etc.
	%

function N = parsenum(N) 
%PARSENUM parse data warehouse mat-file name.
%
%   Example:
%      N = {'out_9_1_10.mat';'out_10_11_1001.mat'};
%      N = repmat(N,100,1);
%

N = strrep(N,'out_',' ');
N = strrep(N,'.mat',';');
N = strrep(N,'_', ',');
N = str2num(cat(2,'[',N{:},']'));
