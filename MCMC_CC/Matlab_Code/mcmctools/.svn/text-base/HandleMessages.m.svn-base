function [] = HandleMessages(varargin)
%HANDLEMESSAGES A error logging program for the MCMC routine
%   This function classifies and logs messages generated during a MCMC
%   search routine. All messages are stored in a global variable LOGS. 
%   There are two ways of calling the function, either with an integer
%   message identifier and a string containing the message, or with a set
%   of string flags. The string can also contain escape sequences and other
%   in line formmating. When sending a message, the identifiers are as 
%   follows:
%       1: Debug message, only should appear if you turn on one or more
%       debug modes. Stored by iteration.
%       2: General message, something which could be relevant for later
%       analysis and interpretation, but doesn't indicate anything is
%       wrong. Stored by iteration.
%       3: Error messages. Store what the error message was that caused the
%       program to crash. Usually only necessary on a distributed or batch
%       job where the program isn't being run on your machine. Stored in
%       one buffer with iteration indicated in message.
%       4: Warning messages. Messages which may indicate something is
%       wrong or unusual, but not necessarily fatal to the routine. Stored
%       in one buffer with iteration indicated in message.
%
%   Each message is automatically timestamped as it is entered in the log.
%   Additionally, one can indicate special variables of interest by adding
%   additional string flags as additonal arguments after the error message,
%   in the form of '-(variable name here)'.
%
%   To call the function using a flag of strings, the following flags may
%   be used:
%   -'iterate': Manually increase the iteration count to one. Note that
%   this program doesn't automatically iterate, to allow multiple messages
%   to be stored in one iteration.
%   -'uniterate': Decrease the iteration count by one. Used in navigating
%   the LOGS variable
%   -'clear': Clears the log for a given iteration.
%
%   

global CONFIG;
global LOGS;
persistent iteration;

now = clock;
now = fix(now);
message = [num2str(now(2)) '/' num2str(now(3)) '/' num2str(now(1)) '-' num2str(now(4)) ':' num2str(now(5)) ':' num2str(now(6)) ' '];
if(nargin >= 2)
    if(mod(varargin{1},1)~= 0|| ~ischar(varargin{2}))
        error('When storing messages, HandleMessages expects an integer identifier and a string to store.');
    end
    if(nargin >2)
        for i = 3:nargin
            if(varargin{i}(1) ~= '-')
                error(['parameter flag ' varargin{i} ' not recognized']);
            else
                message = [message, '(' varargin{i}(2:end) ') '];
            end
        end
    end
    message = [message, varargin{2} , '\n'];
    switch varargin{1}
        case 1
            if(CONFIG.OPTIONS.DEBUG.LogDebug)
                if(isempty(LOGS.Debug{iteration}))
                    LOGS.Debug{iteration} = sprintf(message);
                else
                    LOGS.Debug{iteration} = [LOGS.Debug{iteration}, sprintf(message)];
                end
            end
        case 2
            if(CONFIG.OPTIONS.DEBUG.LogMessages)
                LOGS.Messages{iteration} = [LOGS.Messages{iteration}, sprintf(message)];
            end
        case 3
            if(CONFIG.OPTIONS.DEBUG.LogErrors)
                message = ['ERROR at iteration ' num2str(iteration) ': ' message];
                LOGS.Error = [LOGS.Error, sprintf(message)];
            end
        case 4
            if(CONFIG.OPTIONS.DEBUG.LogWarnings)
                message = ['WARNING at iteration ' num2str(iteration) ': ' message];
                LOGS.Error = [LOGS.Error, sprintf(message)];
            end
        otherwise
            warning(['message identifier ' num2str(varargin{1}) ' not recognized. Message not stored']);
    end
    if(varargin{1} == 1 || varargin{1} == 2)
        disp(sprintf(message));
    end
elseif(nargin == 1)
    if(~ischar(varargin{1}))
        error('When sending one element to HandleMessages, it must be a known flag');
    end
    switch(lower(varargin{1}))
        case 'iterate'
            if(isempty(iteration))
                iteration = 0;
            end
            iteration = iteration + 1;
        case 'uniterate'
            if(iteration > 0)
                iteration = iteration - 1;
            end
        case 'clear'
            LOGS(iteration).Debug = '';
            LOGS(iteration).Messages = '';
        otherwise
            error([varargin{1} ' is an unknown flag.']);
    end
end
            
end

