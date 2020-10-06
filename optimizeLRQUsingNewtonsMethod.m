function [sol_states,u_opt] = optimizeLRQUsingNewtonsMethod(n_rows, n_cols, d, load_mat)
%% Symbolic variables
    x = sym('x', [1, d * n_rows * n_cols], 'real');
    u = sym('u', [1, d * n_rows * n_cols], 'real');
    
    if nargin < 4
        R = createSymbolicMatrix('R', 1, 0, d, 1, 0, 1);
        N = createSymbolicMatrix('N', 1, 0,  d, 1, 0, 1);
        Q = createSymbolicMatrix('Q', 1, 0, d, 1, 0, 1);
        A = createSymbolicMatrix('A', 1, 0, d, 1, 0, 1);
        B = createSymbolicMatrix('B', 1, 0, d, 1, 0, 1);
        C = createSymbolicMatrix('C', 1, 0, d, 1, 0, 1);
        D = createSymbolicMatrix('D', 1, 0, d, 1, 0, 1);
    else
        vars = load(load_mat);
        R = vars.R;
        N = vars.N;
        Q = vars.Q;
        A = vars.A;
        B = vars.B;
        C = vars.C;
        D = vars.D;
    end

    idx_ref = @(i,j) ((j-1)*d+(i-1)*n_cols*d + 1): ((j-1)*d+(i-1)*n_cols*d + d);

%% Cost per state
    cost_per_state = @(i,j) x(idx_ref(i,j)) * Q * x(idx_ref(i,j))' ...
        + u(idx_ref(i,j)) * R * u(idx_ref(i,j))' +...
            2*x(idx_ref(i,j)) * N * u(idx_ref(i,j))';

%% Cost grid
    cost_grid = cell(n_rows,n_cols);
    for i = 1 : n_rows
        for j = 1 : n_cols
            cost_grid{i,j} = cost_per_state(i,j);
        end
    end

%% Assumption grid
    assumption_grid = cell(n_rows,n_cols);
    for i = 1 : n_rows
        for j = 1 : n_cols
            if i == 1 && j == 1
                assumption_grid{i,j} = x(idx_ref(i,j))';
            elseif i == 1
                assumption_grid{i,j} = C * x(idx_ref(i,j-1))' + D * u(idx_ref(i,j))';
            elseif j == 1
                assumption_grid{i,j} = B * x(idx_ref(i-1,j))' + D * u(idx_ref(i,j))';
            else
                assumption_grid{i,j} = A * x(idx_ref(i-1,j-1))' +  D * u(idx_ref(i,j))' + ...
                     B * x(idx_ref(i-1,j))'  + C * x(idx_ref(i,j-1))';
            end
        end
    end
    
    cost = sum([cost_grid{:}]);
    for i = n_rows : -1 : 1
        for j = n_cols : -1 : 1
            cost = subs(cost, x(idx_ref(i,j)), assumption_grid{i,j}');
        end
    end
    
    x_init = ones(1,d);
    cost2 = subs(cost, x(idx_ref(1,1)), x_init);
    initial_u = randn(1, d * n_rows * n_cols);
    
    [grad, hess] = symbolic_gradient_hessian(cost2,u');
    
    [u_opt,num_iter,status_code,status_message] =...
        symbolic_newton_raphson(u', initial_u', grad, hess);
    
    disp(status_message);
    sol_states = cell(n_rows,n_cols);
    
    save('vars_if_failed.mat');
        
    for i = 1 : n_rows
        for j = 1 : n_cols
            if i == 1 && j == 1
                sol_states{i,j} = x_init;
            elseif i == 1 && j > 1
                sol_states{i,j} = double(subs(subs(subs(x(idx_ref(i,j)),...
                    x(idx_ref(i,j)), assumption_grid{i,j}'),...
                    x(idx_ref(i,j-1)),sol_states{i,j-1}), ...
                    u(idx_ref(i,j)),u_opt(idx_ref(i,j))'));
            elseif j==1 && i > 1
                sol_states{i,j} = double(subs(subs(subs(x(idx_ref(i,j)),...
                    x(idx_ref(i,j)), assumption_grid{i,j}'),...
                    x(idx_ref(i-1,j)),sol_states{i-1,j}), ...
                    u(idx_ref(i,j)),u_opt(idx_ref(i,j))'));
            else
                sol_states{i,j} = double(subs(subs(subs(subs(subs(x(idx_ref(i,j)),...
                    x(idx_ref(i,j)), assumption_grid{i,j}'),...
                    x(idx_ref(i-1,j)),sol_states{i-1,j}), ...
                    u(idx_ref(i,j)),u_opt(idx_ref(i,j))'),...
                    x(idx_ref(i,j-1)),sol_states{i,j-1}), ...
                    x(idx_ref(i-1,j-1)),sol_states{i-1,j-1}));
            end
        end
    end
    fprintf('Optimal cost using Newton method is: %f\n', double(subs(cost2, u', u_opt)));
end
