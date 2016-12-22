function ARGS = amb_parse_arguments(argv, options, defaults)

% ARGS = amb_parse_arguments(varargin, default_options)
% ARGS = amb_parse_arguments(varargin, options, defaults)
%
% ARGS - struct
% default_options - cell array of alternating string/value pairs
% options - cell array of strings 
% defaults - cell array of default values (in string form)
% varargin - cell array of function arguments
%
% Returns a struct with fields given by 'options' containing either the
% corresponding values in 'defaults' or the value following a char-array of
% the same option name in varargin, i.e. varargin should be a cell array of
% option string/value pairs set alternately. If the option name appears
% more than once in varargin, the last associated value is used. When the
% two argument form is used, 'default_options' should take the form of
% varargin.

FUNCTION_NAME = 'amb_parse_arguments';

num_args = length(argv);

if nargin==2
  default_options = options;
  num_default_options = length(default_options);
  if mod(num_default_options,1)==1,
    % error is handled below
    num_options = (num_default_options+1)/2;
    num_defaults = num_options - 1;
  else
    num_options = num_default_options/2;
    num_defaults = num_options;
  end
  
  options = {};
  for i = 1:num_options
    options{i} = default_options{2*i-1};
  end
  for i = 1:num_defaults
    defaults{i} = default_options{2*i};
  end
else
  num_options = length(options);
  num_defaults = length(defaults);
end

if num_defaults<num_options
  disp([FUNCTION_NAME ': not enough default values.']);
  for i=num_defaults+1:num_options
    defaults{i} = [];
  end
end

for i=1:num_options
  if ischar(options{i})
    eval(['ARGS.' options{i} ' = defaults{i};']);
  else
    disp([FUNCTION_NAME ': non-char option name requested as option number ' num2str(i)])
  end
end

i = 1;
while i<num_args
  option = argv{i};
  if ischar(option)
    i = i+1;
    if isfield(ARGS,option)
      eval(['ARGS.' option ' = argv{i};']);
    else
      disp([FUNCTION_NAME ': invalid option ''' option ''': ignoring.']);
    end
  else
    disp([FUNCTION_NAME ': expecting char array, found ' num2str(option) ' in argument list: ignoring.'])
  end
  i = i+1;
end



