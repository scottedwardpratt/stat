function OUT = showtr_obj(X,LOGBF,t,FILTER,FLAG,varargin)
%SHOWTR_OBJ Objective function for SHOWTR.
%

persistent the_line PSTRUCT FIELDNAMES TRACE_IGNORE;

if ~isempty(FILTER)
	X = feval(FILTER,X,t,varargin{:});
end


if isstruct(X)
	X = [X(:)]; %Ensure column depth orientation...
	if isempty(PSTRUCT)
		PSTRUCT = logical(1);
		FIELDNAMES = fieldnames(X);
		
		temp = struct2cell(X);
		TRACE_IGNORE = logical(zeros(1,size(temp,1)));
		
		for i = 1:size(temp,1)
			if ~isnumeric(temp{i})
				warning('Trace can not be constructed with non-numeric parameters.')
				TRACE_IGNORE(i) = 1;
				disp(' ');
				disp(['Ignoring ''',FIELDNAMES{i},''' in trace.']);
				disp(' ');
			elseif ~isequal(size(temp{i}),[1 1])
				warning('Trace can not be constructed with vector-valued parameters.')
				TRACE_IGNORE(i) = 1;
				disp(' ');
				disp(['Ignoring ''',FIELDNAMES{i},''' in trace.']);
				disp(' ');
			end
		end

		if ~any(TRACE_IGNORE) %Just to simplify checking later on...
			TRACE_IGNORE = [];
		end
	end
	
	
	if isempty(the_line)
		if ~isempty(TRACE_IGNORE)

			TEMP = struct2cell(X);
			TEMP = TEMP(~TRACE_IGNORE,:);
			TEMP = cell2mat(TEMP)';
			the_line = line(t,TEMP);
			%legend(FIELDNAMES{~TRACE_IGNORE},-1);
			
		else
			TEMP = struct2cell(X);
			TEMP = cell2mat(TEMP)';

			the_line = line(t,TEMP);
			%legend(FIELDNAMES{:},-1);
		end
	else
		
		if ~isempty(TRACE_IGNORE)
			TEMP = struct2cell(X);
			TEMP = TEMP(~TRACE_IGNORE,:);
			TEMP = cell2mat(TEMP)';
			

		else
			TEMP = struct2cell(X);
			TEMP = cell2mat(TEMP)';

		end
		
		for i = 1:length(the_line)
			old_x = get(the_line(i),'XData');
			new_x = [old_x,t];
			old_y = get(the_line(i),'YData');
			new_y = [old_y,TEMP(:,i)'];
			set(the_line(i),'XData',new_x,'YData',new_y);
		end

	end
	
	
else

	if isempty(the_line)
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

if strcmp(FLAG,'report')
	rndlegend(the_line,FIELDNAMES{:});
end

OUT = the_line;

