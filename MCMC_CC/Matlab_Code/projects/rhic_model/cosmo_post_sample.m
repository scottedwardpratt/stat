function [] = cosmo_post( )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global CONFIG;

Results = ProcessResults();

F2 = fieldnames(Results);

for i = 1:size(F2,1)
    disp(['Parameter: ' F2{i} ' mean: ' num2str(Results.(F2{i})(1)) ' standard deviation: ' num2str(Results.(F2{i})(2))]);
end

save('FOOBAR');

end

