THETA = [0.15 0.86 0.83 3.2 1.65 0.48 0.19];

[meanvals, errors] = QueryEmulator('$SCRIPT_HOME/../model-data/rhic-test-run', THETA);

data = rand(1,45); %parse in actual data.

SIGMA = diag(errors);

like = mvnpdf(data-meanvals, 0, SIGMA);

%add in hierarchical errors.
