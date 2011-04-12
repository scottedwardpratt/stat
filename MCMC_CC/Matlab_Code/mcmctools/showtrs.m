function the_line = showtrs(OPTS,varargin)
%SHOWTRS Plot a trace of the progress of MCMC with subplots.
%   SHOWTRS(OPTS) plots traces from the parameter output files
%   in OPTS.OutDir.
%   SHOWTRS(OPTS,N) plots traces using iterations N and greater.
%
%   SHOWTRS(OPTS,I) plots traces using just the iterations in the
%   integer vector I.
%
%   SHOWTRS(OPTS,I,FILT) applies the FILTER prior to plotting the trace.
%
%   H = SHOWTRS(...) returns in H handles to the line objects.
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

the_line = postfun(OPTS,I,'showtrs_obj',FILTER,PARAMS{:});
grid on;
box on;
set(gca,'NextPlot',next_plot_spec);

if ~isempty(FILTER)
	clear(FILTER);
end
	
clear('showtrs_obj');
