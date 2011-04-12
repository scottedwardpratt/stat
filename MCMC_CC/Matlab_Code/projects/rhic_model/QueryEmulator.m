function [ meanvalues, errors ] = QueryEmulator( model_data_file, THETA )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(~exist(model_data_file, 'file'))
    error([model_data_file ' does not exist.']);
end


filename = [getenv('SCRIPT_HOME') '/madaiPts.txt'];

if(isstruct(THETA))
    dlmwrite(filename, (cell2mat(struct2cell(THETA)))','delimiter', '\t', 'precision', '%f');
else
    dlmwrite(filename, THETA ,'delimiter', '\t', 'precision', '%f');
end

CMD = ['cd $SCRIPT_HOME; cat madaiPts.txt | ./computePoints.sh ' model_data_file];

[status,result] = system(CMD);

% 
% if(status ~=0)
%     if(HM), HandleMessages(3, 'Error in running emulator.'); end;
%     error('Error in running emulator.');
% end
% 
% 
% if(DEBUG)
%     if(HM)
%         HandleMessages(1, 'Done.');
%     else
%         disp('Done.');
%     end
% end

foo = zeros(2,45);

if(status == 0)
    s1 = regexp(result, '\n', 'split');
    k=1;
    for i=1:length(s1)
        if(strncmp(s1(i),'#',1) ~= 1 && strcmp(s1(i), '') ~=1)
            if(strncmp(s1(i), 'Error',5) ==1)
                result
                error(s1(i))
            end
            s2 = regexp(s1(i), ' ', 'split');
            s2 = s2{1};
            s3 = str2double(s2);
            
            if(any(isnan(s3))==1)
                result
                s2
                s3
                error('NaNs detected in output.');
            end
            
            try
                foo(k,:) = s3(:);
                k = k+1;
            catch exception
                s2
                s3
                exception
                error('These lines are screwed up.');
            end
        end
    end
else
    error('could not execute emulator.');
end

if (mod(size(foo,1),2)~= 0)
    size(foo,1)
    error('Mismatched number of values and errors. Check emulator output.');
end

meanvalues = zeros(1, 45);
errors = zeros(1, 45);

for i=1:1
    meanvalues(i, :) = foo(2*i-1, :);
    errors(i, :) = foo(2*i, :);
end

delete(filename);
end

