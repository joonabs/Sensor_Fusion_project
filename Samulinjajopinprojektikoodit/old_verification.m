function [theta, thetas, g, R] = verification(y_data, p0, G_symbolic, th_estimate, z, methods)    
    N = size(y_data,2);                          % Number of measurements

    % Measurement model, function handle g(theta)
    g = @(p) kron(ones(N,1), [eye(3) (3*[p; z]*[p; z]'-norm([p; z]).^2*eye(3))/norm([p; z]).^5]*th_estimate);
    
    R = kron(eye(N), cov(y_data'));                 % Measurement noise covariance matrix    
    Imax = 20;                                      % Maximum number of iterations (optional, default: 10)    
                                                    
    y_data = y_data(:);                             % Make y to be 3*Nx1 vector


    % Jacobian matrix of the measurement model, function handle                 
    G_jacobian = @(p) jacobian_substitute(G_symbolic, p, th_estimate(4:6), N);

    [theta, thetas, Js, gammas] = lsqsolve(y_data, g, G_jacobian, R, p0, Imax, method);
end


