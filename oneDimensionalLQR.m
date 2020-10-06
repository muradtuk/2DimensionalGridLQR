function [sol_states, optimal_u_per_sub_array] = oneDimensionalLQR(x, u, R, N, Q, A, B, assumption_array)
n = length(x);
%% Cost per state
    cost_per_state = @(i) x{i}' * Q * x{i} ...
        + u{i}' * R * u{i} +...
            2*x{i}' * N * u{i};

%% Cost array
    cost_array = cell(n,1);
    for i = 1 : n
        cost_array{i} = cost_per_state(i);
    end
 
%% Substitute each x by the assumption and apply the gradient as well
    sub_cost_array = cost_array;
    optimal_u_per_sub_array = cell(n,1);
    for i = n : -1 : 1
        if i == 1 
            break;
        elseif i < n
            sub_cost_array{i,1} = sub_cost_array{i+1,1} + cost_array{i,1};
            sub_cost_array{i,1} = subs(sub_cost_array{i,1}, x{i},...
                assumption_array{i});
        end
        temp = computeOptimalU(sub_cost_array{i,1}, u{i});
        if isstruct(temp)
            optimal_u_per_sub_array{i,1} = cell2sym(struct2cell(...
                temp));
        else
            optimal_u_per_sub_array{i,1} = temp;                
        end
        sub_cost_array{i} = subs(sub_cost_array{i,1}, u{i},...
            optimal_u_per_sub_array{i,1});
    end
    
    %% Fill in the array
    sol_states = x;
    sol_states{1} = subs(sol_states{1}, x{1}, ones(size(x{1})));
    
    for i = 2 : n
        sol_states{i} = subs(subs(assumption_array{i},u{i},...
            optimal_u_per_sub_array{i,1}), x{i-1},...
            sol_states{i-1});
    end    
end

