% Returns a symbolic jacobian of the dynamic model
function G = jacobian_moving_symbolic()
    syms p_x p_y p_z m_x m_y m_z v alfa real
    p = [p_x; p_y; p_z]; p_trans = [p_x p_y p_z];
    m = [m_x; m_y; m_z];
    
    C = @(alfa) [cos(alfa) -sin(alfa) 0;
           sin(alfa) cos(alfa) 0; 0 0 1];
    g = ((3*p*p_trans - sqrt(p_x^2 + p_y^2 + p_z^2)^2*eye(3))/sqrt(p_x^2 + p_y^2 + p_z^2)^5)*C(alfa)*m;
    
    G = jacobian(g, [p_x, p_y, v, alfa]);
end