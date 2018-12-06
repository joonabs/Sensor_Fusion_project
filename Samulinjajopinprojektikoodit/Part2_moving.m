clear all
clc

z = 11;
[th_estimate R] = calibration();
R = 1*eye(3);
m = th_estimate(4:6);

%Y = load_data('Moving_circle_data_1.txt');       % Measurement data
Y = load_data('Moving_circle_data_2.txt');       % Measurement data
N = size(Y,2);

%x0 = [30 0 1 3/2*pi]';
x0 = [-25 5 1 3/2*pi]';
P0 = diag([sqrt(5) sqrt(5) 0.1 0.5]);


sigma_v = 0.01;
sigma_alfa = 0.1;
t_delta = 0.1;

Q = Q_symbolic(sigma_v, sigma_alfa, t_delta);
F_symbolic = quasi_jacobian_symbolic(t_delta);
G_symbolic = jacobian_moving_symbolic();

G = @(x) jacobian_moving_substitute(G_symbolic, x, m);
F = @(x) quasi_jacobian_substitute(F_symbolic, x);

C = @(alfa) [cos(alfa) -sin(alfa) 0; sin(alfa) cos(alfa) 0; 0 0 1];
g = @(p, alfa) [eye(3) (3*[p; z]*[p; z]'-norm([p; z]).^2*eye(3))/norm([p; z]).^5]* ...
        [eye(3) zeros(3); zeros(3) C(alfa)]*th_estimate;

f = @(x) x + t_delta*x(3)*[cos(x(4)) sin(x(4)) 0 0]';
    
    
x = x0; P = P0;
estimates = []; P_values = []; K_values = [];
for i = 1:N
    % Prediction
    %F = quasi_jacobian_substitute(F_symbolic, x, t_delta*i);
    x = f(x);
    P = F(x)*P*F(x)' + Q_substitute(Q, x);
    
    Gx = G(x);
    
    % Measurement update
    K = P*Gx'/(Gx*P*Gx' + R);
    x = x + K*(Y(:,i) - g(x(1:2), x(4)));
    P = P - K*(Gx*P*Gx' + R)*K';
    
    P_values = [P_values P];
    K_values = [K_values K];
    estimates = [estimates x];
end

figure()
hold on
plot(x0(1), x0(2), 'xr')
plot(estimates(1,:), estimates(2,:), 'b')
legend('Starting point','Estimated')
axis equal
