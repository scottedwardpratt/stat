function OUT = showtrs_obj(X,t,FILTER,FLAG,varargin)
%SHOWTR_OBJ Objective function for SHOWTRS.
%

persistent the_ax the_line;

if ~isempty(FILTER)
	X = feval(FILTER,X,t,varargin{:});
end


if (isempty(the_line) | isempty(the_ax))
	row = floor(sqrt(size(X,2)));
	col = ceil(size(X,2)/row); 
	
	for i = 1:size(X,2)	
		the_ax(i) = subplot(row,col,i);
		the_line(i) = line(t,X(:,i));	
	end 

else

	for j = 1:length(the_line)
		subplot(the_ax(i)); %Make ith subplot current.
		old_x = get(the_line(j),'XData');
		new_x = [old_x,t];
		old_y = get(the_line(j),'YData');
		new_y = [old_y,X(:,j)'];
		set(the_line(j),'XData',new_x,'YData',new_y);
	end

end

OUT = the_line;

