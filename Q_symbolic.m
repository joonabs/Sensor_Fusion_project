function Q = Q_symbolic(siqma_v, siqma_alfa, t_delta)
    syms tn tn_previous tau p_x p_y v alfa w_x w_y real

    I = eye(4);
    B = [0 0; 0 0; 1 0; 0 1];
    w = [w_x; w_y];
    siqma = [siqma_v 0; 0 siqma_alfa];

    f = [v*cos(alfa); v*sin(alfa); 0; 0] + B*w;
    A = jacobian(f, [p_x, p_y, v, alfa]);

    q = (I + A*(tn - tau))*B*siqma*B'*(I + A*(tn - tau))';

    Q = int(q, tau, [tn_previous tn]);
    tn_previous = tn - t_delta;
    Q = subs(Q);  
end