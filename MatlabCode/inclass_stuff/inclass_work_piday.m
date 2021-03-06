% Pi Day lecture notes
%% Cell arrays and structures, Chp 4.1

a1 = ['Joe', 'Schmoe', 'Netflix'] % This horizontally concats them, output = JoeSchmoeNetflix

% a2 = ['Joe'; 'Schmoe'; 'Netflix';] % This throws an error, must be equal size to vertically stack

a3 = ['Joe    '; 'Schmoe '; 'Netflix'] % This properly vertically stacks them

% Cell arrays dont need to worry about size

a4 = {'Joe', 'Schmoe', 'Netflix'} % Returns a 1x3 cell array

a5 = {'Joe'; 'Schmoe'; 'Netflix';} % Returns a 3x1 cell array

e = {2, 'steve'; 1:2, 1,2; 3,4}; % Returns a 2x2 cell array

e{2,2}(3) % This returns the number 2

e{1,2} % Returns steve

celldisp(e) % displays every part of the cell in an odd graphic

%% Functions of variable length, 4.2
clear, clc, close all

s = perimeter('square', 5)
r = perimeter('rectangle', 4, 3)
c = perimeter('circle', 10)

% varargout function calls
% Perimeter2 is the new version
Ps1 = perimeter2('square', 5)
Pr1 = perimeter2('rectangle', 4, 3)
Pc1 = perimeter2('circle', 10)
[Ps, As]= perimeter2('square', 5)
[Pr, Ar] = perimeter2('rectangle', 4, 3)
[Pc, Ac] = perimeter2('circle', 10)

function output = perimeter(varargin) 
    % Putting varargin will do different things based on the amount of arguments passed to it (variable)
    % Preset the output to 0
    output = 0;
    
    % Make a switch statement to check the first input (which will be 'square', 'circle', etc
    switch varargin{1}
        case 'square'
            % Check the number of inputs given ('square',5) has nargin of 2
            if nargin == 2 
                % varargin is a cell array of inputs, so our length is at {2}
                length = varargin{2};
                output = 4*length; 
            end
            
        case 'rectangle'
            if nargin == 3
                % two lengths input, so {2} and {3} are l, w
                length = varargin{2};
                width = varargin{3};
                output = 2*length + 2*width; 
            end
            
        case 'circle'
            % cell entry {2} is our radius
            r = varargin{2};
            output = 2*pi*r;
    end
end

function varargout = perimeter2(varargin) 
    % Putting varargin will do different things based on the amount of arguments passed to it (variable)
    % Preset the output to 0
    varargout{1} = 0;
    varargout{2} = 0;
    
    % Make a switch statement to check the first input (which will be 'square', 'circle', etc
    switch varargin{1}
        case 'square'
            % Check the number of inputs given ('square',5) has nargin of 2
            if nargin == 2 
                % varargin is a cell array of inputs, so our length is at {2}
                length = varargin{2};
                varargout{1} = 4*length; 
                varargout{2} = length*length;
            end
            
        case 'rectangle'
            if nargin == 3
                % two lengths input, so {2} and {3} are l, w
                length = varargin{2};
                width = varargin{3};
                varargout{1} = 2*length + 2*width; 
                varargout{2} = length*width;
            end
            
        case 'circle'
            % cell entry {2} is our radius
            r = varargin{2};
            varargout{1} = 2*pi*r;
            varargout{2} = pi*r*r;
    end
end