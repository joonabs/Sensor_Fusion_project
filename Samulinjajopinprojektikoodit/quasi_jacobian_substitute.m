% Substitutes values of p and m to jacobian and returns numerical answer
function G = quasi_jacobian_substitute(G, x)
    v = x(3);
    alfa = x(4);
    %t = current_t;

    G = subs(G);
    G = double(G);
end