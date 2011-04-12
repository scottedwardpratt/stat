function [ string ] = struct2comment( s, type )
%struct2comment converts a structure s into a comment, whose style is
%specified by variable type. 
%   
%   Suppose we have a complex structure like the CONFIG variable generated
%   by Configure.m.
%   Example:
%   
%   global CONFIG;
%   Configure('set');
%   struct2comment(CONFIG, 'latex');
%

if(~ischar(type))
    error('type must be a character string');
end

name = inputname(1);

%Show what type is being displayed for reimporting back into Matlab


switch type
    case 'text'
        string  = sprintf(['#DISPTYPE = ' type '\n']);
        %Display the variable name
        string = [string, sprintf(['#VARIABLE = ' name '\n'])];
        %get the structure's internals with seperate function-call
        structureguts = structcrawl(s,type,{name});
        string = [string, structureguts];
        
    case 'latex'
        string  = sprintf('%s\n',['\#DISPTYPE = ' type ' \\']);
        PREAMBLE = '';
        PREAMBLE = [PREAMBLE,sprintf('%s\n','\documentclass{article}')];
        PREAMBLE = [PREAMBLE,sprintf('%s\n','\usepackage{amsmath}')];
        PREAMBLE = [PREAMBLE,sprintf('%s\n','\usepackage{pstricks}')];
        PREAMBLE = [PREAMBLE,sprintf('%s\n','\begin{document}')];
        PREAMBLE = [PREAMBLE,sprintf('%s\n','\setlength{\parindent}{0em}')];
        

        string = [string, sprintf('%s\n', ['\#VARIABLE = ' name])];
        string = [string, sprintf('%s\n', '\begin{itemize}')];
        %add slashes to underscores so latex can process them
        name = strrep(name, '_','\_');
        %in case the name was already latex formatted and we added a second
        %slash
        name = strrep(name, '\\_','\_');
        structureguts = structcrawl(s,type,{name});
        string = [string, structureguts];

        string = [PREAMBLE,string];

        DOCEND = [];
        DOCEND = [DOCEND,sprintf('%s\n','\end{itemize}')];
        DOCEND = [DOCEND,sprintf('%s\n','\end{document}')];

        string = [string,DOCEND];

    otherwise
        error([type ' is an unknown display type.']);
end

end

function  result = structcrawl( s, type, parent )
%A recursive funtion that performs a depth first search to display the
%names and values of all the structure's variables.

F = fieldnames(s);
Vals = struct2cell(s);
result = '';
for i = 1:length(F)    
    switch type
        case 'text'
            %Determine what type of data the value is, and translate
            %appropriately
            Class= class(Vals{i});
            switch Class
                %T/F: stored as 'true' and 'false' as opposed to Matlab's
                %1/0 convention.
                case 'logical'
                    if(Vals{i})
                        result = [result, sprintf(['#TYPE ' Class ' ' F{i} ' = true\n'])];
                    else
                        result = [result, sprintf(['#TYPE ' Class ' ' F{i} ' = false\n'])];
                    end
                %cell arrays: depending on what the cell array contains,
                %translate differently. NOTE: known bug: this assumes that
                %the cell array contains all one type of data, and will
                %error if it encounters another type.
                case 'cell'
                    if(~isempty(Vals{i}))
                        switch class(Vals{i}{1})
                            case 'double'
                                %for a cell array of values, concatenate into a
                                %normal array, convert to string, and print
                                result = [result, sprintf(['#TYPE cell double ' F{i} ' = ' num2str([Vals{i}{:}]) '\n'])];
                            case 'char'
                                %Can't concatenate a character cell array,
                                %spaces are lost. Therefore, print one at a
                                %time with spaces inserted
                                tempstr = ['#TYPE cell char ' F{i} ' = '];
                                for j = 1:length(Vals{i})
                                    tempstr = [tempstr, Vals{i}{j} ' '];
                                end
                                tempstr = [tempstr, '\n'];
                                result = [result, sprintf(tempstr)];
                            otherwise
                                error([F{i} ' contains a cell array of unknown type ' class(Vals{i}{1})]);
                        end
                    else
                        result = [result, sprintf(['#TYPE cell empty ' F{i} ' = '])];
                    end
                %Numbers: convert to string and print    
                case 'double'
                    result = [result, sprintf(['#TYPE double ' F{i} ' = ' num2str(Vals{i}) '\n' ])];
                %Structures: Display #CATEGORY tag for reimporting, then
                %recusively call structcrawl again.
                case 'struct'
                    result = [result, sprintf(['#CATEGORY ' F{i} ' of '])];
                    for j  = 1 : length(parent)
                        result = [result, sprintf([parent{j} ' '])];
                    end
                    result = [result,sprintf('\n')];
                    parents = [parent ,F(i)];
                    result = [result, structcrawl(Vals{i},type,parents)];
                    result = [result, sprintf('#END \n')];
                case 'char'
                    %In filepaths, the '\' characters must be converted to
                    %'\\' in order to avoid erroneous escape sequences.
                    tempstr = strrep(Vals{i}, '\' , '\\');
                    
                    result = [result, sprintf(['#TYPE char ' F{i} ' = ' tempstr ' \n'])];
                otherwise
                    error([F{i} ' has unknown type ' Class]);
            end
        case 'latex'
            %Replace underscores with escape clauses so latex can process
            %them, and correct for the extra slash if F{i} was already
            %latex formatted.
            F{i} = strrep(F{i}, '_', '\_');
            F{i} = strrep(F{i}, '\\_', '\_');
            
            Class= class(Vals{i});
            
            switch Class
                %T/F: stored as 'true' and 'false' as opposed to Matlab's
                %1/0 convention.
                case 'logical'
                    if(Vals{i})
                        result = [result, sprintf('%s\n',['\item \#TYPE ' Class ' ' F{i} ' = true'])];
                    else
                        result = [result, sprintf('%s\n',['\item \#TYPE ' Class ' ' F{i} ' = false'])];
                    end
                %cell arrays: depending on what the cell array contains,
                %translate differently. NOTE: known bug: this assumes that
                %the cell array contains all one type of data, and will
                %error if it encounters another type.
                case 'cell'
                    if(~isempty(Vals{i}))
                        switch class(Vals{i}{1})
                            case 'double'
                                %for a cell array of values, concatenate into a
                                %normal array, convert to string, and print
                                result = [result, sprintf('%s\n', ['\item \#TYPE cell double ' F{i} ' = ' num2str([Vals{i}{:}])])];
                            case 'char'
                                %Can't concatenate a character cell array,
                                %spaces are lost. Therefore, print one at a
                                %time with spaces inserted
                                tempstr = ['\item \#TYPE cell char ' F{i} ' = '];
                                for j = 1:length(Vals{i})
                                    Vals{i}{j} = strrep(Vals{i}{j}, '_','\_');
                                    Vals{i}{j} = strrep(Vals{i}{j}, '\\_', '\_');
                                    tempstr = [tempstr, Vals{i}{j} ' '];
                                end
                                result = [result, sprintf('%s\n',tempstr)];
                            otherwise
                                error([F{i} ' contains a cell array of unknown type ' class(Vals{i}{1})]);
                        end
                    else
                        result = [result, sprintf('%s\n', ['\item \#TYPE cell empty ' F{i} ' = '])];
                    end
                %Numbers: convert to string and print    
                case 'double'
                    result = [result, sprintf('%s\n',['\item \#TYPE double ' F{i} ' = ' num2str(Vals{i})])];
                %Structures: Display #CATEGORY tag for reimporting, then
                %recusively call structcrawl again.
                case 'struct'
                    temp = sprintf('%s', ['\item \#CATEGORY ' F{i} ' of ']);
                    for j  = 1 : length(parent)
                        temp = [temp, sprintf([parent{j} ' '])];
                    end
                    result = [result,sprintf('%s\n',temp)];
                    result = [result, sprintf('%s\n', '\begin{itemize}')];
                    parents = [parent ,F(i)];
                    result = [result, structcrawl(Vals{i},type,parents)];
                    result = [result, sprintf('%s\n','\end{itemize}')];
                case 'char'
                    %In filepaths, the '\' characters must be converted to
                    %'\\' in order to avoid erroneous escape sequences.
                    tempstr = strrep(Vals{i}, '\' , '\\');
                    tempstr = strrep(tempstr, '_', '\_');
                    tempstr = strrep(tempstr, '\\_', '\_');
                    
                    result = [result, sprintf('%s\n',['\item \#TYPE char ' F{i} ' = ' tempstr])];
                otherwise
                    error([F{i} ' has unknown type ' Class]);
            end
            
        otherwise
            error([type ' is an unknown display type.']);
    end
end

end
