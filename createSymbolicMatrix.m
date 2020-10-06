function variable = createSymbolicMatrix(var_name, n_rows, n_cols,d, is_symmetric, is_vector, is_not_cell)
    variable = cell(n_rows,n_cols);
    for i = 1 : n_rows
        for j = 1 : n_cols
            temp_var_name = strcat(var_name, num2str(i), '_',...
                    num2str(j), '_');
            if is_vector
                variable{i,j} = sym(temp_var_name, [d 1], 'real');
            else
                if is_symmetric
                    temp = triu(sym(temp_var_name, [d, d], 'real'));
                    temp = temp + temp' - diag(diag(temp));
                    variable{i,j} = temp;
                else
                    variable{i,j} = sym(temp_var_name, [d, d], 'real');
                end
            end
        end
    end
    
    if n_rows == 1 && is_not_cell == 1
        variable = randi(9, d,d);
    end
end



function variable = createSymbolicMatrixOLD(var_name, n,d, is_symmetric, is_vector)
    
    if is_vector
        variable = sym(var_name, [n n d], 'real');
    else
        variable = sym(var_name, [n n d d], 'real');
        if is_symmetric
            for i = 1 : n
                for j = 1 : n
                    temp = triu(reshape(variable(i,j, :, :), [d, d]));
                    temp = temp + temp' - diag(diag(temp));
                    variable(i,j,:,:) = temp;
                end
            end        
        end
    end
end

