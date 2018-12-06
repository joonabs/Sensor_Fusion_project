% Substitutes values of "WHICH ONES" to jacobian and returns a numerical answer
function G = jacobian_moving_substitute(G, x, m)
    p_x = x(1);
    p_y = x(2);
    p_z = 11;
    v = x(3);
    alfa = x(4);
    
    m_x = m(1);
    m_y = m(2);
    m_z = m(3);

    G = subs(subs(G));
    G = double(G);
end