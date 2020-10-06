function optimal_u = computeOptimalU(cost_func,variable)
    grad = gradient(cost_func, variable);
    eqn = (grad == 0);
    optimal_u = solve(eqn, variable);
end

