function variable = createSymbolicArray(var_name, n,d, is_symmetric, is_vector, multiple_paths)
    if nargin < 6
        multiple_paths = 1;
    end
    variable = cell(multiple_paths,n);
    for path = 1 : multiple_paths
        for i = 1 : n
            temp_var_name = strcat(var_name,'_', num2str(path), '_', num2str(i));
            if is_vector
                variable{path, i} = sym(temp_var_name, [d 1], 'real');
            else
                if is_symmetric
                    temp = triu(sym(temp_var_name, [d, d], 'real'));
                    temp = temp + temp' - diag(diag(temp));
                    variable{path,i} = temp;
                else
                    variable{path,i} = sym(temp_var_name, [d, d], 'real');
                end
            end
        end 
    end
    
    if n == 1
        if multiple_paths == 1
            variable = randi(9, d,d);
        else
            for i = 1 : multiple_paths
                variable{i} = randi(9, d,d);
            end
        end
    end
end

