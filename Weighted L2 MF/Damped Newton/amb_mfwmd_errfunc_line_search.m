function [x2,t2,fx2,exit_code] = amb_mfwmd_errfunc_line_search(x1,dir,AUX,varargin)

% [x_new,alpha,error,success] = 
%   amb_mfwmd_errfunc_line_search(x_start,direction,AUX,options...)
%
% Performs a line search in 'direction' from 'x_start' until a minimum is 
% found, where upon that new x is returned. 'alpha' is such that "x_new = 
% x_start + alpha * direction". 'error' is the error evaluated at 'x_new'.
%
% 'AUX' must be a struct with the following fields (all with usual
% meanings): m, n, r, M and W 
%
% 'success' is 1 if a minimum is found, zero otherwise (if the line search
% iteration limit was reached).
%
% Options:
%
% 'input_checks'   [1]   check inputs are consistent
% 'sumsqrderr'     [1]   coeff for amb_mfwmd_errfunc_summsqrderr
% 'regularizeAB'   [0]   coeff for amb_mfwmd_errfunc_regularizeAB 
% 'regratio'       [1]   ratio of reg(B) to reg(A) in ...regularizeAB 
% 'line_search_tolerance' [1e-4] line search termination threshold
% 'line_search_limit' [20]  maximum number of line search steps


FUNCTION_NAME = 'amb_mfwmd_errfunc_line_search';

default_options = { ...
    'input_checks', 1, ...  
    'sumsqrderr', 1, ...
    'regularizeAB', 0, ...
    'regratio', 1, ...
    'line_search_tolerance', 1e-4, ...
    'line_search_limit', 20, ...
  };


ARGS = amb_parse_arguments(varargin,default_options);

if ARGS.input_checks
  % default return values
  x2 = x1;
  t2 = 0;
  fx2 = nan;
  exit_code = 1;
  
  if ~(isfield(AUX,'m') && isfield(AUX,'n') && isfield(AUX,'r') && isfield(AUX,'M') && isfield(AUX,'W'))
    disp([FUNCTION_NAME ': AUX must have the fields m, n, r, M and W.']);
    return
  end
  
  if ~(isequal(size(AUX.M),[AUX.m,AUX.n]) && isequal(size(AUX.W),[AUX.m,AUX.n]))
    disp([FUNCTION_NAME ': AUX.M and AUX.W must be m-by-n.']);
    return
  end
  
  if length(x1)~=AUX.r*(AUX.m+AUX.n) || length(dir)~=AUX.r*(AUX.m+AUX.n)
    disp([FUNCTION_NAME ': x_start and direction must be r*(m+n)-vectors.']);  
    return
  end
end
  
fx_initial = ambfwmderrfunc_vec(x1,AUX,ARGS);
fx1 = fx_initial;

% parameterize on multiples of the gradient vector
t_step_factor = 1.2;
t_step = 1;
t1 = 0;
t2 = t1 + t_step; t_step = t_step*t_step_factor;
t3 = t2 + t_step; t_step = t_step*t_step_factor;

% find other side of local minimum (minima)
x2 = x1 + t2*dir; fx2 = ambfwmderrfunc_vec(x2,AUX,ARGS);
x3 = x1 + t3*dir; fx3 = ambfwmderrfunc_vec(x3,AUX,ARGS);

% if fx2 is higher than fx1 then go and do a bisection step, otherwise
% slide along line until fx3 is higher than fx2.
if fx2<fx1
  while fx3<fx2
    t1 = t2; t2 = t3; t3 = t2 + t_step;
    
    x1 = x2; fx1 = fx2;
    x2 = x3; fx2 = fx3;
    x3 = x2 + t_step*dir; fx3 = ambfwmderrfunc_vec(x3,AUX,ARGS);
    
    t_step = t_step*t_step_factor;
  end  
end

% quadratic bisection line search
% (x1 and x3 straddle (hopefully only one) minimum)
line_search_improvement = inf; error_improvement=inf; line_search_count = 1;
while error_improvement>ARGS.line_search_tolerance & line_search_improvement>ARGS.line_search_tolerance & line_search_count<=ARGS.line_search_limit
  last_fx2 = fx2;
  
  % create matrix N with rows (t^2,t,1,-y) to fit quadratic. However, N
  % is likely to be ill-conditioned as max(abs([fx1 fx2 fx3]))>>t3 is
  % very probable. Fortunately, scaling will not alter the quadratic's 
  % stationary point. 
  scale = max(abs([fx1 fx2 fx3]));
  N = [t1^2 t1 1 -fx1/scale; t2^2 t2 1 -fx2/scale; t3^2 t3 1 -fx3/scale];
  q = null(N); % a = q(1)/q(4); b = q(2)/q(4); c = q(3)/q(4);
  t = -0.5*q(2)/q(1); % t = -b/(2*a);
  
  if q(1)/q(4)>=0 % quadratic with minimum \_/
    line_search_improvement = abs(t-t2);
  
    if t<t1
      % t<--t1 t2 t3
      x1 = x1+(t-t1)*dir; fx1 = ambfwmderrfunc_vec(x1,AUX,ARGS); t1 = t;
      last_fx2 = inf; % fx2 not actually moved
    elseif t>t3 
      % t1 t2 t3-->t
      x3 = x3+(t-t3)*dir; fx3 = ambfwmderrfunc_vec(x3,AUX,ARGS); t3 = t;
      last_fx2 = inf; % fx2 not actually moved
    elseif t-t2>ARGS.line_search_tolerance
      % t1-->t2-->t t3
      t1 = t2; fx1 = fx2; x1 = x2;
      x2 = x2+(t-t2)*dir; fx2 = ambfwmderrfunc_vec(x2,AUX,ARGS); t2 = t;
    elseif t-t2<-ARGS.line_search_tolerance
      % t1 t<--t2<--t3
      t3 = t2; fx3 = fx2; x3 = x2;
      x2 = x2+(t-t2)*dir; fx2 = ambfwmderrfunc_vec(x2,AUX,ARGS); t2 = t;
    end % else tolerance achieved
  
    line_search_count = line_search_count + 1;
  else % quadratic with maximum /^\ : slide
    t_step = (t3-t2)*t_step_factor;
    
    while fx3<fx2
      t1 = t2; t2 = t3; t3 = t2 + t_step;
      
      x1 = x2; fx1 = fx2;
      x2 = x3; fx2 = fx3;
      x3 = x2 + t_step*dir; fx3 = ambfwmderrfunc_vec(x3,AUX,ARGS);
      
      t_step = t_step*t_step_factor;
    end
    line_search_count = line_search_count + 0.01;
  end  
  error_improvement = abs(last_fx2 - fx2); % hmmm - maybe increase in fx2 should be noticed...
end

if fx2>fx_initial, fx2 = fx_initial; end % final failsafe

exit_code = line_search_count<=ARGS.line_search_limit;

function e = ambfwmderrfunc_vec(x,AUX,ARGS)
  A = reshape(x(1:AUX.m*AUX.r),AUX.r,AUX.m)';
  B = reshape(x(AUX.m*AUX.r+1:end),AUX.r,AUX.n)';
  e = 0;
  if ARGS.sumsqrderr>0, e = e + ARGS.sumsqrderr*amb_mfwmd_errfunc_sumsqrderr(A,B,AUX.M,AUX.W); end
  if ARGS.regularizeAB>0, e = e + ARGS.regularizeAB*amb_mfwmd_errfunc_regularizeAB(A,B,ARGS.regratio); end
return

