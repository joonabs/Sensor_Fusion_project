%   Projectwork dynamic
%
%
%   Aleksi Mäkinen and Joonas Isometsä
clear
clc

%   Initial parameters
%   
%   Height of the magnet
z = 10;

%   Initialization parameters of spectral density of white noise processes
sigma_v = 0.01;
sigma_alfa = 0.1;
t_delta = 0.1;

%   Initial guess for the starting point of movement
%   And the corresponding known actual positions:
%   [x1, y1, z1; x2, y2, z2; ... etc]
p = [30,0,10; 15,25,10; 0,45,10; -30,45,10; -30,0,10; -20,-30,10; 15,-25,10; 0,-55,10]';

%   Circular movement 1:
%data_filu_nimi = 'circular_movement_1.txt';
%x0 = [30; 0; 10; pi/2];
%P0 = diag([1 1 0.1 0.5);

%   Straightline movement:
data_filu_nimi = 'J_A_straight_line_movement.txt';
x0 = [43, -33, 10, pi/2]';
P0 = diag([1 1 0.1 0.5]);

%   First calibrate the model
%   y               contains measurements used for calibration (used just to see the
%                   situation
%   th_estimate     contains the estimated superposition of dipole moment and
%                   earths magnetic field
[y, th_estimate] = calibration(p);

% The estimate is now renamed m_T !
m_T = th_estimate;

%   Open measurement data of the movement
y_data = load_data(data_filu_nimi);
%   Save the amount of measurement points
N = size(y_data, 2);

%   Calculate and save symbolic Q, F and G
Q = Q_symbolic(sigma_v, sigma_alfa, t_delta);

F = @(x) [
    1, 0, t_delta*cos(x(4)), -t_delta*x(3)*sin(x(4));
    0, 1, t_delta*sin(x(4)),  t_delta*x(3)*cos(x(4));
    0, 0,            1,                  0;
    0, 0,            0,                  1;
];

G_symbolic = jacobian_symbolic_G();

C = @(alfa) [cos(alfa) -sin(alfa) 0;...
    sin(alfa) cos(alfa) 0;...
    0 0 1];
g = @(p, alfa) [eye(3) (3*[p; 10]*[p; 10]' - norm([p; 10]).^2*eye(3))/...
    norm([p; 10]).^5]*[eye(3) zeros(3); zeros(3) C(alfa)]*m_T;

% Input initial values
x = x0;
P = P0;
% Preallocate room for estimates
x_estimatelist = zeros(4, N);
P_estimatelist = zeros(4, 4, N);

% The main loop:
for i = 1:N
    % Prediction
    % Discretized quasi-constant model estimates next step
    x = x + [t_delta*x(3)*cos(x(4));...
        t_delta*x(3)*sin(x(4)); 0; 0];
    
    P = F(x)*P*F(x)' + Q_substitute(Q, x);    
    
    G = jacobian_substitute(G, x, m_T);
    
    % Measurement update
    K = P*G'/(G*P*G' + R);
    x = x + K*(y_data(:,i) - g(x(1:2), x(4)));
    P = P - K*(G*P*G' + R)*K';
    
    x_estimatelist(:, i) = x;
    P_estimatelist(:, i) = P;
end

% Plots of the results

figure()
hold on
plot(estimates(1,:), estimates(2,:))
legend('Estimated')


figure()
hold on
%plot(X(1,:), 'b-')
plot(estimates(1,:), 'r-')
title('x-coordinate')
legend('Estimated')

figure()
hold on
%plot(X(2,:), 'b-')
plot(estimates(2,:), 'r-')
title('y-coordinate')
legend('Estimated')



%{
%   Solved symbolic F Jacobian
function F = jacobian_symbolic_F()
    syms alfa v
    F = [1 0 cos(alfa) -v*sin(alfa);...
        0 1 sin(alfa) v*cos(alfa);...
        0 0 1 0;...
        0 0 0 1];
end
%}
%   Symbolic G jacobian:
function G = jacobian_symbolic_G()
    syms p_x p_y p_z m_T_1 m_T_2 m_T_3 alfa real
    
    p = [p_x ; p_y ; p_z];
    m_T = [m_T_1; m_T_2; m_T_3];
    C = [cos(alfa) -sin(alfa) 0;...
        sin(alfa) cos(alfa) 0;...
        0 0 1];
    
    y = ((3*p*p' - norm(p).^2*eye(3))/norm(p).^5)*C*m_T
    G = jacobian(y, [p_x, p_y]);
    
end
%   Function to insert and calculate numeric jacobian
function G = jacobian_substitute(G, x, m_T)
    p_x = x(1);
    p_y = x(2);
    p_z = 10;
    m_T_1 = m_T(4);
    m_T_2 = m_T(5);
    m_T_3 = m_T(6);
    alfa = x(4);
    
    G = subs(G);
    G = double(G);
end

function Q = Q_substitute(G, x)
    % tn tn_previous tau p_x p_y v alfa
    p_x = x(1);
    p_y = x(2);
    p_z = 10;
    alfa = x(4);
end
