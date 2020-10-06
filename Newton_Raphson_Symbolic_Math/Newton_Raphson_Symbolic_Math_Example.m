% Example: Newton Raphson Optimization (Minimization) using symbolic math
% demonstrates the functions 'symbolic_newton_raphson' and 'symbolic_gradient_hessian'
% Please refer to the functions for complete documentation.
% 
% Written by Dr. Yoash Levron, Technion, Israel, 2015

clc;

%%%%% step I - create a symbolic function %%%%

% Choose an example
select_example = 1;    % See 'switch' statement below
N = 10;  % select number of variables

X = sym('X',[N 1]);   % create symbolic vector of N variables
p = (1:N)'; % parameters of the symbolic function
% in the following examples f is a symbolic function which global minimum is X=p
switch select_example
    case 1
        % This function is quadratic in X, so the algorithm will converge to
        % the global optimum with a single iteration
        f = sum((p+1).*(X-p).^2);
    case 2
        % This function is convex, so the algorithm will converge to
        % the global optimum after a few iterations
        f = sum((X-p).^4);
    case 3
        % The Hessian matrix of this matrix is singular, and an infinite
        % number of (global) minimum points exists. The algorithm converges
        % to an arbitrary minimum point at a single iteration
        f = sum(diag(X-p)*ones(N,N)*(X-p));        
    case 4
        % This function is not convex so the Newton algorithm may converge to a local minimum
        f = sum(  (cos(X-p)+(X-p)-1).^2  );
    otherwise
        disp('Unknown example');
        return
end

%%%%% Step II - Compute the gradient vector & Hessian matrix:
[grad hess] = symbolic_gradient_hessian(f,X);

%%%%% Step III- Find the minimum point of 'f' by the Newton Raphson Algorithm.
Xinit = 0.5*ones(N,1); % guess an arbitrary initial solution
% Call the Newton Raphson solver:
[xopt,num_iter,status_code,status_message] = symbolic_newton_raphson(X,Xinit,grad,hess);

% Display results
clc
disp(status_message);
optimal_solution_X = xopt
number_of_iterations = num_iter
function_value_at_solution  = double(subs(f,X,xopt))

