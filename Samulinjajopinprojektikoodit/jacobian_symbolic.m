% Calculates Jacobian symbolically
function G = jacobian_symbolic()
    syms p_x p_y p_z
    syms m_x m_y m_z
    p = [p_x; p_y; p_z];
    p_trans = [p_x p_y p_z];
    m = [m_x; m_y; m_z];
    
    g = ((3*p*p_trans - sqrt(p_x^2 + p_y^2 + p_z^2)^2*eye(3))/sqrt(p_x^2 + p_y^2 + p_z^2)^5)*m;
    G = jacobian(g, [p_x, p_y]);
end