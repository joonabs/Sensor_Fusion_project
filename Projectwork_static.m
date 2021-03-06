% Project work - static part
% Aleksi M�kinen and Joonas Isomets�
clear all
clc
%   The initial parameters required for the code to work:
%      
%   For the code to work input:
%   Measurements to calibrate the model from measurements-file, at least 2 needed!
%   In the form: 'sensorLog_%i.txt'

%   The corresponding known actual positions:
%   [x1, y1, z1; x2, y2, z2; ... etc]
% vanhat pisteet p = [30, 0, 10; 30, 30, 10; 0, 45, 10; 30, -30, 10]';
% alla uudet pisteet
p = [30, 0, 10; 
    15, 25, 10;
    10, 45, 10;
    -30, 45, 10;
    -30, 0, 10;
    -20, -30, 10 ;
    15, -25, 10]';
%   Measurements to validate the model from measurements-file, at least 2 needed!
%file21 = 'sensorLog_20181025T140033.txt';
file21 = 'J_A_static_4.txt';
%file21 = 'sensorLog_20181025T133655.txt';
%   And the actual position (so that the accuracy can be compared in
%   plot!):
%actual_pos = [15, 30, 10]';
% pointti 8 ei ollenkaan k�yt�ss� ja estimoidaan pistett� 1
actual_pos = [-30, 45 10]';
%actual_pos = [45,15,10]';
%   Initial position guess for validation:
%   [x0, y0, z0]'
theta0 = [-10, 10]';
%   Maximum number of iterations of the LSQsolve (optional, def 10):
Imax = 40;
%    Method to use, one of 'gradient', 'gauss-newton', or
%    'levenberg-marquardt' (optional, default: 'gradient')
%method = 'gradient';
%method = 'gauss-newton';
method = 'levenberg-marquardt';

%   First the calibration
%
%   Use the given load_data function:
%   ym      3xNm matrix of magnetometer measurement data. The rows are the
%           three magnetometer axes (x, y, z), the columns the samples in 
%           time. 
%           The rest aren't used now.
%[ym1, tm, ya, ta] = load_data(file11);

%   Assemble the y and the respective G:s
y = [];
G = [];
for i = 1:size(p,2)
    filename = sprintf('J_A_static_%i.txt',i);
    data = load_data(filename);
    y_current = data(:);
    
    p_current = p(:,i);
    H = (3*p_current*p_current'-norm(p_current).^2*eye(3))/norm(p_current).^5; 
    G_current = [eye(3) H];

    G = [G; repmat(G_current, size(y_current,1) / 3, 1)];
    y = [y; y_current];
end

% Now the estimate can be calculated:
th_estimate = (G'*G)\G'*y;

%       Then the validation
%
%       The following values and functions are needed:
%
%   y       Measurement data (3N x 1)
%   g       Measurement model, function handle g(theta)
%   G       Jacobian matrix of the measurement model, function handle 
%           G(theta)
%   R       Measurement noise covariance matrix (3 x 3)
%   theta0  Already set
%   Imax    - || -
%   method  - || -
%       
%   We start by loading the validation data (3 x N):
[ym1, tm, ya, ta] = load_data(file21);

%   Save the amount of measurements 3*N
N = size(ym1, 2);


%   Now the measurement noise covariance matrix R (3 x 3):
R = cov(ym1');

%   Scale it to suitable form for LSQsolve (3N x 3N):
R = kron(eye(N), R);

%   The measurement modelfunction handle (scaled from :
g = @(p) kron(ones(N,1), [eye(3) (3*[p; 10]*[p; 10]' - norm([p; 10]).^2*eye(3))/norm([p; 10]).^5]*th_estimate);

    %   Change ym1 to suitable form for LSQsolve (3N x 1):
ym = ym1(:);

%   Jacobian function handle
% aleksin vanha ratkasu Jacobian = @(p) G_jacobian_new(p, N, th_estimate(4:6));
Jacobian = jacobian_symbolic();
Jacobian =  @(p) jacobian_substitute(Jacobian, p, th_estimate, N);


%   Finally it's time to estimate the location using given function:

[theta, thetas, Js, gammas] = lsqsolve(ym, g, Jacobian, R, theta0, Imax, method);

%   Next we can plot the initial guess and the results
%   First the cost function for the plotting:
J = @(p) (ym - g(p))'/R*(ym - g(p)); 

%   Next construct the grid for the plot:
[pxx, pyy] = meshgrid(-50:1:0, 0:1:50);

Js = zeros(size(pxx));
for i = 1:size(pxx, 1)
    for j = 1:size(pyy, 2)
        Js(i, j) = J([pxx(i, j); pyy(i, j)]);
    end
end

figure(1); clf();
contourf(pxx, pyy, log(Js), 40);                % from example
%contour(th1grid, th2grid, costFunc, 100);  % from exercise 3.1

hold on;

%   Plot the actual end position
plot(actual_pos(1), actual_pos(2), 'g*');
%   Plot final result
plot(theta(1), theta(2), 'or');
%   Plot the advancement of the algorithm
plot(thetas(1, :), thetas(2, :), '-xr');
%   Plot the initial guess
plot(theta0(1), theta0(2), '+k')

%axis([-60 75, -60, 75]);
axis equal;

legend('Cost function','Target true location', 'Final result', 'Progressing line' ,'Initial guess')

if isequal(method, 'gradient')
    title('Gradient Descent')
end 

if isequal(method, 'levenberg-marquardt')
    title('Levenberg-Marquardt')
end 

if isequal(method, 'gauss-newton')
    title('Gauss-Newton')
end 


function G = jacobian_symbolic()
    syms p_x p_y p_z m_T_1 m_T_2 m_T_3 real;
    
    p = [p_x ; p_y ; p_z];
    m_T = [m_T_1; m_T_2; m_T_3];
    y = ((3*p*p' - norm(p).^2*eye(3))/norm(p).^5)*m_T;
    G = jacobian(y, [p_x, p_y]);
    
end

function G = jacobian_substitute(G, x, th_estimate, N)
    p_x = x(1);
    p_y = x(2);
    p_z = 10;
    m_T_1 = th_estimate(4);
    m_T_2 = th_estimate(5);
    m_T_3 = th_estimate(6);
    
    
    G = subs(G);
    G = double(G);
    G_stack = [];
    for i = 1:N
        G_stack = [G_stack ; G];
    end 
    G = G_stack;
end

%   Testing purposes:
%G_theta = Jacobian(theta0)
