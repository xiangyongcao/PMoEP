function [A,B,states] = amb_mfwmd(A,B,M,W,varargin)

% [A,B,states] = amb_mfwmd(A,B,M,W,options...)
%
% Matrix factorization with missing data.
%
% A range of algorithms to minimize 
% e(A,B) = norm(W.*(M-A*B'))^2 + L * ( norm(A)^2 + Lr*norm(B)^2 )
%
% Options:
%
% 'algorithm'       ['damped Newton'] algorithm to use for minimization
% 'tolerance'       [1e-10]   stop when error decreases by less than tol.
% 'iterations'      [1e3]   maximum number of iterations
% 'L'               [1e-7]  regularizing magnitude
% 'Lr'              [1]     regularizer ratio
% 'lambda'          [1]     initial lambda in damped Newton
% 'lambda_step'     [8]     the shift applied in damped Newton
% 'micro_tolerance' [1e-6]  "brief" state records error at this improvement
% 'brief'           [0]     set to record final error and number iterations
% 'errors'          [0]     set to one to record total errors
% 'separate'        [0]     set to one to record both error components
% 'times'           [0]     set to one to record timings
% 'matrices'        [0]     set to one to record all As and Bs 
% 'switches'        [0]     set to one to record hybrid alt/dn switches
% 'graph'           [0]     set to one to have real-time graphing
%
% Algorithms
%
% 'aanaes'                  'Robust Factorization' Aanaes et al, PAMI:24(9)
% 'altdn01'                 alternation/damped newton hybrid version 01
% 'altdn01a'                alternation/damped newton hybrid version 01a
% 'altdn02'                 alternation/damped newton hybrid version 02
% 'altdn02a'                alternation/damped newton hybrid version 02a
% 'altdn03'                 alternation/damped newton hybrid version 03
% 'alternation'             the basic alternation algorithm
% 'brandt'                  'Closed-Form ...' Brandt ECCV'02 Workshop
% 'damped Newton'           the damped Newton algorithm
% 'dn altreg'               the 'altreg' variation of damped Newton
% 'dn line search'          damped Newton with line search
% 'dn ls altreg             'altreg' damped Newton with line search
% 'gradient descent'        standard gradient descent
% 'huynh'                   'Outlier Correction ...' Huynh et al. ICCV'03
% 'powerfactorization'      'PowerFactorization' Vidal and Hartley CVPR'04
% 'project and merge'       the obvious rank truncation scheme
% 'shum'                    'PCA with missing data...' Shum et al PAMI:19(9)
%
% states is a struct with some or all of the following fields: 
%   algorithm - name and function name of the algorithm used
%   total_errors - arrary of the total error for each iteration
%   sumsqrd_errors - array holding just norm(W.*(M-A*B'))^2
%   regularizing_errors - array of only L*(norm(A)^2 + norm(B)^2)
%   times - array with the cumulative time for each iteration
%   matrices - struct array holding the As and Bs every iterations
%   switches - for which iterations alternation was used (hybrids only)
%   final_error - error after last iteration
%   iterations - total number of iterations
%   micro_error - the error after improvement less than 'micro_tolerance'
%   micro_iteration - the iteration when micro_error recorded
%   final_improvement - the difference in error of last two iterations

FUNCTION_NAME = 'amb_mfwmd';

options = {'algorithm','tolerance','iterations','L','Lr','lambda','lambda_step','micro_tolerance','brief','errors','separate','times','matrices','switches','graph'};
defaults = {'damped Newton',1e-10,1e2,1e-7,1,0.1,8,1e-6,0,0,0,0,0,0,0};

ARGS = amb_parse_arguments(varargin,options,defaults);

states = struct;

if ARGS.errors || ARGS.graph, states.total_errors = nan*ones(ARGS.iterations,1); end
if ARGS.separate, states.sumsqrd_errors = nan*ones(ARGS.iterations,1); end
if ARGS.separate, states.regularizing_errors = nan*ones(ARGS.iterations,1); end
if ARGS.times, states.times = nan*ones(ARGS.iterations,1); end
if ARGS.matrices, states.matrices = struct('A',[],'B',[]); end

if ARGS.brief
    states.tolerance = ARGS.tolerance;
    states.micro_tolerance = ARGS.micro_tolerance;
    states.micro_error = nan;
    states.micro_error_iteration = nan;
    micro_error_recorded = 0;
end

if ARGS.graph
    FUNCTION_NAME_SPECIAL = [FUNCTION_NAME ' (' ARGS.algorithm ')']; 
    FUNCTION_NAME_SPECIAL(find(FUNCTION_NAME_SPECIAL=='_')) = ' ';
    figure(gcf)
    % line_h = semilogy(nan,'k');
end

[m,n,r,valid] = amb_mfwmd_matrix_dimension_check(A,B,M,W,FUNCTION_NAME);
if ~valid, return, end

k = 1;
last_error = inf;
this_error = amb_mfwmd_errfunc_sumsqrderr(A,B,M,W) + ARGS.L * amb_mfwmd_errfunc_regularizeAB(A,B,ARGS.Lr);

% alg_state initialization
alg_state = struct('A',A,'B',B,'M',M,'W',W,'tolerance',ARGS.tolerance,'L',ARGS.L,'Lr',ARGS.Lr,'lambda',ARGS.lambda,'lambda_step',ARGS.lambda_step,'error',this_error);

switch lower(ARGS.algorithm)
    case 'aanaes'
        algorithm_fhandle = @amb_mfwmd_aanaes;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'altdn01'
        algorithm_fhandle = @amb_mfwmd_alt_dn_hybrid_v01;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r,ARGS.switches); % extra alg_state initialization
    case 'altdn01a'
        algorithm_fhandle = @amb_mfwmd_alt_dn_hybrid_v01a;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r,ARGS.switches); % extra alg_state initialization
    case 'altdn02'
        algorithm_fhandle = @amb_mfwmd_alt_dn_hybrid_v02;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r,ARGS.switches); % extra alg_state initialization
    case 'altdn02a'
        algorithm_fhandle = @amb_mfwmd_alt_dn_hybrid_v02a;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r,ARGS.switches); % extra alg_state initialization
    case 'altdn03'
        algorithm_fhandle = @amb_mfwmd_alt_dn_hybrid_v03;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r,ARGS.switches); % extra alg_state initialization
    case 'alternation'
        algorithm_fhandle = @amb_mfwmd_alternation;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'brandt'
        algorithm_fhandle = @amb_mfwmd_aanaes;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'damped newton'
        algorithm_fhandle = @amb_mfwmd_damped_newton;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'dn altreg'
        algorithm_fhandle = @amb_mfwmd_dn_altreg;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'dn line search'
        algorithm_fhandle = @amb_mfwmd_dn_line_search;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'dn ls altreg'
        algorithm_fhandle = @amb_mfwmd_dn_ls_altreg;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'gradient descent'
        algorithm_fhandle = @amb_mfwmd_gradient_descent;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'huynh'
        algorithm_fhandle = @amb_mfwmd_huynh;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'powerfactorization'
        algorithm_fhandle = @amb_mfwmd_powerfactorization;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    case 'project and merge'
        algorithm_fhandle = @amb_mfwmd_project_and_merge;
        alg_state = feval(algorithm_fhandle,alg_state, r); % extra alg_state initialization
    case 'shum'
        algorithm_fhandle = @amb_mfwmd_shum;
        alg_state = feval(algorithm_fhandle,alg_state, m,n,r); % extra alg_state initialization
    otherwise
        %disp([FUNCTION_NAME ': unknown algorithm [' ARGS.algorithm '] specified.'])
        error(['Specified algorithm [' ARGS.algorithm '] is unknown.'])
        algorithm_fhandle = 0;
        ARGS.iterations = 0; % quit immediately, but output state correct
end

if isa(algorithm_fhandle, 'function_handle'), states.algorithm = [ARGS.algorithm ' (' func2str(algorithm_fhandle) ')'];
else states.algorithm = [ARGS.algorithm ' (?)'];
end

if ARGS.errors || ARGS.graph, states.total_errors(k) = this_error; end
if ARGS.separate, states.sumsqrd_errors(k) = amb_fwmd_errfunc_sumsqrderr(A,B,M,W); end
if ARGS.separate, states.regularizing_errors(k) = ARGS.L * amb_fwmd_errfunc_regularizeAB(A,B,ARGS.Lr); end
if ARGS.times, states.times(k) = 0; end
if ARGS.matrices, states.matrices(k).A = A; states.matrices(k).B = B; end

improvement = last_error - this_error;

tic;

k = k + 1;
while k<=ARGS.iterations && improvement>ARGS.tolerance
	last_error = this_error;
	
	% iterate desired minimization algorithm
	alg_state = feval(algorithm_fhandle,alg_state);
	
	this_error = alg_state.error;
    improvement = last_error - this_error;
	
	if ARGS.errors || ARGS.graph, states.total_errors(k) = this_error; end
	if ARGS.separate, states.sumsqrd_errors(k) = amb_fwmd_errfunc_sumsqrderr(alg_state.A,alg_state.B,M,W); end
	if ARGS.separate, states.regularizing_errors(k) = ARGS.L * amb_fwmd_errfunc_regularizeAB(alg_state.A,alg_state.B,ARGS.Lr); end
	if ARGS.times, states.times(k) = states.times(k-1)+toc; tic, end
	if ARGS.matrices, states.matrices(k).A = A; states.matrices(k).B = B; end
    if ARGS.brief && not(micro_error_recorded) && improvement<ARGS.micro_tolerance, states.micro_error = this_error; states.micro_error_iteration = k; micro_error_recorded = 1; end
	
	if ARGS.graph
		hold off
		semilogy(states.total_errors(1:k),'k')
        if ARGS.separate
            hold on
            plot(states.sumsqrd_errors(1:k),'r')
            plot(states.regularizing_errors(1:k),'b')
            legend('total','sumsqrd','regularizing')
        end
		title([FUNCTION_NAME_SPECIAL ': iteration ' num2str(k)])
		xlabel(['iteration (improvement = ' num2str(improvement) ')'])
		ylabel('error')
		drawnow
	end
	
	k = k + 1;  
end

almost_instantaneous_time = toc;

A = alg_state.A;
B = alg_state.B;

if ARGS.errors, states.total_errors = states.total_errors(1:k-1); end
if ARGS.separate, states.sumsqrd_errors = states.sumsqrd_errors(1:k-1); end
if ARGS.separate, states.regularizing_errors = states.regularizing_errors(1:k-1); end
if ARGS.times, states.times = states.times(1:k-1); end

if ARGS.graph && not(ARGS.errors), states = rmfield(states,'total_errors'); end

if ARGS.brief
    states.final_error = this_error;
    states.iterations = k-1;
    states.final_improvement = improvement;
end

if ARGS.switches
    if isfield(alg_state, 'switches'), states.switches = alg_state.switches;
    else warning([FUNCTION_NAME ': algorithm did not return switch information.']);
    end
end

