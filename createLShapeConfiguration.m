function variable = createLShapeConfiguration(nH,nV,d)
    variable = cell(n,n);
    for i = 1 : n
        for j = 1 : n
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
    if n == 1
        variable = variable{1};
    end
    
    if n == 1
        variable = randi(9, d,d);
    end
end

