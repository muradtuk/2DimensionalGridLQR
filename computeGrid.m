function [sol_states, optimal_u_per_sub_grid] = computeGrid(n_rows, n_cols, d, load_mat)
%% Symbolic variables
    x = createSymbolicMatrix('x', n_rows, n_cols, d, 0, 1, 0);
    u = createSymbolicMatrix('u', n_rows, n_cols, d, 0, 1, 0);
    
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
    

%% Cost per state
%     cost_per_state = @(i,j) x{i,j}' * x{i,j} ...
%         + u{i,j}' * u{i,j} +...
%             2*x{i,j}' * u{i,j};
    
        
%     cost_per_state = @(i,j) x{i,j}' * Q{i,j} * x{i,j} ...
%         + u{i,j}' * R{i,j} * u{i,j} +...
%             2*x{i,j}' * N{i,j} * u{i,j};
    cost_per_state = @(i,j) x{i,j}' * Q * x{i,j} ...
        + u{i,j}' * R * u{i,j} +...
            2*x{i,j}' * N * u{i,j};

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
                assumption_grid{i,j} = x{i,j};
            elseif i == 1
                assumption_grid{i,j} = C * x{i,j-1} + D * u{i,j};
            elseif j == 1
                assumption_grid{i,j} = B * x{i-1,j} + D * u{i,j};
            else
                assumption_grid{i,j} = A * x{i-1,j-1} +  D * u{i,j} + ...
                     B * x{i-1, j}  + C * x{i, j-1};
            end
        end
    end
 
%% Substitute each x by the assumption and apply the gradient as well
    sub_cost_grid = cell(n_rows,n_cols);
    optimal_u_per_sub_grid = cell(n_rows,n_cols);
    for i = n_rows : -1 : 1
        for j = n_cols : -1 : 1
            if i ==1 && j == 1
                break;
            end
            if j == n_cols && i == n_rows
                sub_cost_grid{i,j} = subs(cost_grid{i,j}, x{i,j}, ...
                    assumption_grid{i,j});
            elseif i == n_rows
                sub_cost_grid{i,j} = ...
                    subs(sub_cost_grid{i,j+1} + cost_grid{i,j}, x{i,j}, ...
                        assumption_grid{i,j});
            elseif j == n_cols
                sub_cost_grid{i,j} = ...
                    subs(sub_cost_grid{i+1,j} + cost_grid{i,j}, x{i,j}, ...
                        assumption_grid{i,j});
            else
                sub_cost_grid{i,j} = ...
                    sub_cost_grid{i+1,j} + sub_cost_grid{i,j+1} + cost_grid{i,j};
                substracted_term = subs(subs(sub_cost_grid{i+1, j+1}, ...
                    x{i, j+1}, assumption_grid{i, j+1}), x{i+1,j}, ...
                        assumption_grid{i+1, j});
                substracted_term = subs(subs(substracted_term, u{i+1,j},...
                    optimal_u_per_sub_grid{i+1,j}), u{i,j+1},...
                    optimal_u_per_sub_grid{i, j+1});
                sub_cost_grid{i,j} = vpa(subs(sub_cost_grid{i,j}...
                    - substracted_term, x{i,j}, assumption_grid{i,j}), 5);
                
            end
            sub_cost_grid{i,j} = clearUnresolvedDependencies(sub_cost_grid{i,j},...
            x,i,j,n_rows);
            temp = computeOptimalU(sub_cost_grid{i,j}, u{i,j});
            if isstruct(temp)
                optimal_u_per_sub_grid{i,j} = vpa(cell2sym(struct2cell(...
                    temp)),5);
            else
                optimal_u_per_sub_grid{i,j} = vpa(temp,5);                
            end
            sub_cost_grid{i,j} = vpa(subs(sub_cost_grid{i,j},...
                    u{i,j}, optimal_u_per_sub_grid{i,j}),5);
        end
    end
    
    %% Fill in the grid
    init_state = ones(size(x{1,1}));
    optimal_u_per_sub_grid{1,1} = zeros(size(u{1,1}));
    sol_states = x;
    
    for i = 1 : n_rows
        for j = 1 : n_cols
            if i == 1
                if j > 1
                    optimal_u_per_sub_grid{i,j} = ...
                        double(subs(optimal_u_per_sub_grid{i,j}, x{i,j-1},...
                        sol_states{i, j-1}));
                    sol_states{i,j} = subs(subs(assumption_grid{i,j},...
                        u{i,j}, optimal_u_per_sub_grid{i,j}), x{i,j-1},...
                        sol_states{i, j-1});
                else
                    sol_states{i,j} = subs(sol_states{i,j}, x{i,j},...
                        init_state); 
                end
            elseif j == 1
                optimal_u_per_sub_grid{i,j} = ...
                        double(subs(optimal_u_per_sub_grid{i,j}, x{i-1,j},...
                        sol_states{i-1, j}));
                sol_states{i,j} = subs(subs(assumption_grid{i,j},...
                        u{i,j}, optimal_u_per_sub_grid{i,j}), x{i-1,j},...
                        sol_states{i-1, j});
            else
                optimal_u_per_sub_grid{i,j} = ...
                        double(subs(subs(subs(optimal_u_per_sub_grid{i,j}, x{i-1,j},...
                        sol_states{i-1, j}), x{i, j-1},...
                        sol_states{i, j-1}), x{i-1,j-1},...
                        sol_states{i-1,j-1}));
                sol_states{i,j} = subs(subs(subs(subs(assumption_grid{i,j},...
                        u{i,j}, optimal_u_per_sub_grid{i,j}), x{i-1,j},...
                        sol_states{i-1, j}), x{i, j-1},...
                        sol_states{i, j-1}), x{i-1,j-1},...
                        sol_states{i-1,j-1});
            end
            sol_states{i,j} = double(subs( sol_states{i,j}, x{1,1},  sol_states{1,1}));
        end
    end
    
    
%% Computing cost
    cost = sum([cost_grid{:}]);
    
    for i = n_rows : -1 : 1
        for j = n_cols : -1 : 1
            cost = subs(subs(cost, u{i,j}, ...
                optimal_u_per_sub_grid{i,j}), x{i,j}, sol_states{i,j}); 
        end
    end
    
    fprintf('Optimal cost is: %f\n', double(cost));    
end



%     sub_cost_grid = cost_grid(n,n);
%     for i = n : -1 : 1
%         for j = n : -1 : 1
%             optimal_u_per_state{i,j} = computeOptimalU(...
%                 sum(sum([sub_cost_grid{i:n, j:n}])), ...
%                     u{i,j});
%         end
%     end


