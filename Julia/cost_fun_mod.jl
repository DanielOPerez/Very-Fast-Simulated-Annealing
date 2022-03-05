module cost_fun_mod
export cost_fun

function cost_fun(var_in::Array,
                  args_in)

    x=var_in[1]
    return sum(100.0*(x[2:end]-x[1:end-1].^2.0).^2.0+(1.0.-x[1:end-1]).^2.0)

end

end
