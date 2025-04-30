
function ini_TaylorN(dim::Int, U::DataType)

    variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:dim+3]
    TaylorSeries.set_variables(U, variable_names, order = 2)
    
end

