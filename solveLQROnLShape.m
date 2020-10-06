function opt_states = lqrRoey(n, d, max_epochs,load_mat)
    %% Symbolic variables
    x = createSymbolicArray('x', n, d, 0, 1, 2);
    u = createSymbolicArray('u', n, d, 0, 1, 2);
    
    if nargin < 4
        R = createSymbolicArray('R', 1, d, 1, 0);
        N = createSymbolicArray('N', 1, d, 1, 0);
        Q = createSymbolicArray('Q', 1, d, 1, 0);
        A = createSymbolicArray('A', 1, d, 1, 0, 2);
        B = createSymbolicArray('B', 1, d, 1, 0, 2);   
    else
%         vars = load(load_mat);
%         R = vars.R;
%         N = vars.N;
%         Q = vars.Q;
%         A = vars.A;
%         B = vars.B;
        load(load_mat);
    end

%% Cost per state
    cost_per_state = @(i,j) x{i,j}' * Q * x{i,j} ...
        + u{i,j}' * R * u{i,j} +...
            2*x{i,j}' * N * u{i,j};

%% Cost grid
    cost_grid = cell(2,n);
    for i = 1 : 2
        for j = 1 : n
            cost_grid{i,j} = cost_per_state(i,j);
        end
    end

%% Assumption grid
    assumption_grid = cell(2,n);
    for path = 1 : 2
        for i = 1 : n
            if i == 1
                assumption_grid{path, i} = x{path,i};
            elseif i == n
                assumption_grid{path, i} = ...
                    A{path} * x{path, i-1} + ...
                    B{path} * u{path, i} + ...
                    A{path + (-1)^(path-1)} *...
                    x{path + (-1)^(path-1), i-1} ...
                    + B{path + (-1)^(path-1)} * u{path + (-1)^(path-1), i};
            else
                assumption_grid{path, i} = A{path} * x{path, i-1}...
                    +  B{path} * u{path,i};
            end
        end
    end
 
%% Substitute each x by the assumption and apply the gradient as well
    opt_states = cell(2,n);
    optimal_u_per_sub_grid = cell(2,n);
    
    for epoch = 1 : max_epochs
        for path = 1 : 2
            if path == 1 && epoch == 1
                temp_assump = assumption_grid(path + 1, 1:end);
                temp_assump{n} = subs(subs(temp_assump{n}, x{path,n-1},...
                    zeros(size(x{path,n-1}))), u{path, n},...
                    zeros(size(u{path,n-1})));
                [opt_states(path+1, 1:end),...
                optimal_u_per_sub_grid(path+1, 1:end)] =...
                oneDimensionalLQR(x(path+1, 1:end),...
                u(path+1, 1:end), R, N, Q, A{path+1},...
                B{path+1}, temp_assump);
            end
            temp_assump = assumption_grid(path, 1:end);
            temp_assump{n} = subs(subs(temp_assump{n},...
                x{path + (-1)^(path-1), n-1}, ...
                opt_states{path + (-1)^(path-1), n-1}),...
                u{path + (-1)^(path-1), n}, ...
                optimal_u_per_sub_grid{path + (-1)^(path-1), n});
            [opt_states(path, 1:end),...
            optimal_u_per_sub_grid(path, 1:end)] =...
            oneDimensionalLQR(x(path, 1:end),...
            u(path, 1:end), R, N, Q, A{path},...
            B{path}, temp_assump);
        end
    end
end

