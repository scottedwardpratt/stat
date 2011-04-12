ACTUAL.x11=NaN;

clear Actual_Plot;
clear Data_Plot;
clear Data_Error;

for i = 1:length(F)
    Actual_Plot{i} = ACTUAL.(F{i});
    Data_Plot{i}=Results.(F{i})(1);
    Data_Error{i}=Results.(F{i})(2);
end
Labels = {'', F{:}, ''};

Actual_Plot = cell2mat(Actual_Plot);
Data_Plot = cell2mat(Data_Plot);
Data_Error = cell2mat(Data_Error);

h1 = plot(Actual_Plot, 'dr');
hold on;
h2 = errorbar(Data_Plot, Data_Error, 'sb');
hold off;
set(gca, 'XTick', 0:1:6);
set(gca, 'XTickLabel', Labels);
title('Comparison of parameter values');
xlabel('Parameter');
ylabel('Parameter value');
legend('Actual Values','Calculated Values');