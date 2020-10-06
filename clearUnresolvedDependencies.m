function refined_cost = clearUnresolvedDependencies(cost,x,i,j,n)
    refined_cost = cost;
    
    if i > 1
        for i_prime = i-1 : -1 : 1
            for j_prime = j+1 : n
                refined_cost = subs(refined_cost, x{i_prime,j_prime},...
                    zeros(size(x{i_prime,j_prime}))); 
            end
        end
    end
    
    if j > 1
        for i_prime = i+1 : n
            for j_prime = j-1 : -1 : 1
                refined_cost = subs(refined_cost, x{i_prime,j_prime},...
                    zeros(size(x{i_prime,j_prime}))); 
            end
        end
    end
end

