function the_line = showtr(OPTS,varargin)
%SHOWTR Plot a trace of the progress of MCMC.
%   SHOWTR(OPTS) plots a trace from the parameter output files
%   in OPTS.OutDir.
%   SHOWTR(OPTS,N) plots a trace using iterations N and greater.
%
%   SHOWTR(OPTS,I) plots a trace using just the iterations in the
%   integer vector I.
%
%   SHOWTR(OPTS,I,FILT) applies the FILTER prior to plotting the trace.
%
%   H = SHOWTR(...) returns in H handles to the line objects.
%

%Keep in mind that it will typically not be possible to hold all the
%data in memory at once.
%
I = [];
if (nargin >= 2)
	I = varargin{1};
end

FILTER = [];
PARAMS = {};

if (nargin >= 3)
	FILTER = varargin{2};
	PARAMS = {varargin{3:end}};
end

if isempty(OPTS.OutDir)
	error('OutDir is empty.  No data can be found.');
end

next_plot_spec = get(gca,'NextPlot');
newplot;


the_line = postfun(OPTS,I,'showtr_obj',FILTER,PARAMS{:});
grid on;
box on;
set(gca,'NextPlot',next_plot_spec);

if ~isempty(FILTER)
	clear(FILTER);
end
	
clear('showtr_obj');

% S = dir(OPTS.OutDir);
% 
% if (length(S) == 0)
% 	error('Can''t get guesses -- no OUT files in OutDir.');
% end
% 
% N = strvcat(S.name);
% ind_out = find(strcmp(cellstr(N(:,1:4)),'out_'));
% %Get the numbers...
% [N,ind] = sort(str2num(N(ind_out,5:end-4)));
% ind_out = ind_out(ind);
% 
% t = 0;
% 
% next_plot_spec = get(gca,'NextPlot');
% newplot;
% hold on;
% 	
% for i = 1:length(ind_out)
% 	X = getfield(load(fullfile(OPTS.OutDir,S(ind_out(i)).name),'X'),'X');
% 	t = [t(end)+1:t(end)+size(X,1)];
% 	
% 	
% 	if (i == 1)
% 		the_line = line(t,X);
% 		
% 	else
% 		for j = 1:length(the_line)
% 			old_x = get(the_line(j),'XData');
% 			new_x = [old_x,t];
% 			old_y = get(the_line(j),'YData');
% 			new_y = [old_y,X(:,j)'];
% 			set(the_line(j),'XData',new_x,'YData',new_y);
% 		end
% 	end
% 
% end
% 
% set(gca,'NextPlot',next_plot_spec);
