% Calculates Quasi jacobian symbolically
function G = quasi_jacobian_symbolic(t_delta)
    %syms p_x p_y v alfa t real
    syms p_x p_y v alfa real
    x = [p_x p_y v alfa]';
    
    f = x + t_delta*x(3)*[cos(x(4)) sin(x(4)) 0 0]';
    
    G = jacobian(f, [p_x, p_y, v, alfa]);
end