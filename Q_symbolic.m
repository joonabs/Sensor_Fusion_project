function Q = Q_symbolic(siqma_v, sigma_alfa)
    syms tn tn_previous tau p_x p_y v alfa w_x w_y real

    I = eye(4);
    B = [0 0; 0 0; 1 0; 0 1];
    w = [w_x; w_y];
    sigma = [sigma_v 0; 0 sigma_alfa];

    f = [v*cos(alfa); v*sin(alfa); 0; 0] + B*w;
    A = jacobian(f, [p_x, p_y, v, alfa]);

    q = (I + A*(tn - tau))*B*sigma*B'*(I + A*(tn - tau))';

    Q = int(q, tau, [tn_previous tn]);
end