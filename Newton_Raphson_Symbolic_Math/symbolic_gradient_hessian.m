function [grad hess] = symbolic_gradient_hessian(f,X)
% Computes the gradient vector & Hessian matrix of a symbolic function
% Inputs:
% f - symbolic function f(X)
% X - symbolic vector of variables
% Outputs
% grad - symbolic gradient vector
% hess - symbolic Hessian matrix
%
% Written by Dr. Yoash Levron, Technion, Israel, 2015

% Gradient vector:
N = size(X,1);
grad = sym('temp',[N 1]);   % initialize
for ii=1:N
    grad(ii) = diff(f,X(ii));
    grad(ii) = simplify(grad(ii));
end
% Hessian matrix:
hess = sym('temp',[N N]);
for ii=1:N
    disp('computing derivatives ...');
    disp([ii N]);
    for jj=1:N
        hess(ii,jj) = diff(grad(ii),X(jj));
        hess(ii,jj) = simplify(hess(ii,jj));
    end
end
disp('done computing derivatives');

end

