function Q = Q_substitute(Q, x)
    v = x(3);
    alfa = x(4);

    Q = subs(Q);
    Q = double(Q);
end