function H = showtr(OPTS)
%SHOWTR Plot a trace of the progress of MCMC.
%   SHOWTR(OPTS) plots a trace from the parameter output files
%   in OPTS.OutDir.
%

%Keep in mind that it will typically not be possible to hold all the
%data in memory at once.
%


S = dir(OPTS.OutDir);

if (length(S) == 0)
	error('Can''t get guesses -- no OUT files in OutDir.');
end

N = strvcat(S.name);
ind_out = find(strcmp(cellstr(N(:,1:4)),'out_'));
%Get the numbers...
[N,ind] = sort(str2num(N(ind_out,5:end-4)));
ind_out = ind_out(ind);

t = 0;

next_plot_spec = get(gca,'NextPlot');
newplot;
hold on;
	
for i = 1:length(ind_out)
	X = getfield(load(fullfile(OPTS.OutDir,S(ind_out(i)).name),'X'),'X');
	t = [t(end)+1:t(end)+size(X,1)];
	
	
	if (i == 1)
		the_line = line(t,X);
		
	else
		for j = 1:length(the_line)
			old_x = get(the_line(j),'XData');
			new_x = [old_x,t];
			old_y = get(the_line(j),'YData');
			new_y = [old_y,X(:,j)'];
			set(the_line(j),'XData',new_x,'YData',new_y);
		end
	end

end

set(gca,'NextPlot',next_plot_spec);
