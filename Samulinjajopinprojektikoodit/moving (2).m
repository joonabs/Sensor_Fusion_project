clear
clc

x0 = [0 0 1 pi/2]';
P0 = eye(4);
t = 0.01;                       % 100 Hz meaning 100 measurements per second

th_estimate = calibration();

C = @(alfa) [cos(alfa) -sin(alfa) 0; sin(alfa) cos(alfa) 0; 0 0 1];
g = @(p, alfa) kron(ones(N,1), [eye(3) (3*[p; z]*[p; z]'-norm([p; z]).^2*eye(3))/norm([p; z]).^5]*C(alfa)*th_estimate);

F_jacobian = quasi_jacobian_symbolic(t);
F = @(x) quasi_jacobian_substitute(F_jacobian, x);

G_jacobian = jacobian_moving_symbolic()
G = @(x) jacobian_moving_substitute(G_jacobian, x, th_estimate);


Q = eye(4);
R = eye(4);

X = [];

f = @(x) x + t*x(3)*[cos(x(4)) sin(x(4)) 0 0]';
x = x0;
P = P0;

for i = 1:1000
    % Prediction
    x = f(x);
    X(:,i) = x;                 % Save model locations
    P = F(x)*P*F(x)' + Q;

    % Measurement
    K = P*G(x)'/(G(x)*P*G(x)' + R);
    x = x + K*(y(:,i) - g(x));
    P = P - K*(G(x)*P*G(x)' + R)*K';
end

figure();
plot(X(1,:),X(2,:))
axis equal