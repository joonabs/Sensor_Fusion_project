% Substitutes values of p and m to jacobian and returns numerical answer
function G = jacobian_substitute(G, p, m, N)
    p_x = p(1);
    p_y = p(2);
    %p_z = p(3);
    p_z = 11;
    
    m_x = m(1);
    m_y = m(2);
    m_z = m(3);
    
    G = subs(G);
    G = double(G);
    G = kron(ones(N,1), G);
end