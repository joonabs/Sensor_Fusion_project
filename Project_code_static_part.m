% Project work - static part
% Aleksi Mäkinen and Joonas Isometsä
clear
clc
%   The initial parameters required for the code to work:
%      
%   For the code to work input:
%   Measurements to calibrate the model from measurements-file, at least 2 needed!
%   In the form: 'sensorLog_%i.txt'

%   The corresponding known actual positions:
%   [x1, y1, z1; x2, y2, z2; ... etc]
p = [30, 0, 10; 30, 30, 10; 0, 45, 10]';

%   Measurements to validate the model from measurements-file, at least 2 needed!
file21 = 'sensorLog_4.txt';

%   Initial position guess for validation:
%   [x0, y0, z0]'
theta0 = [25,-25,10]';
%   Maximum number of iterations of the LSQsolve (optional, def 10):
Imax = 100;
%    Method to use, one of 'gradient', 'gauss-newton', or
%    'levenberg-marquardt' (optional, default: 'gradient')
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
    filename = sprintf('sensorLog_%i.txt',i);
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
[k1, k2] = size (ym1);
N = k2;

%   Change it to suitable form for LSQsolve (3N x 1):
ym = [];
for i = 1 : k2
    current_y = ym1(:,i);
    ym = [ym; current_y];
end

%   Now the measurement noise covariance matrix R (3 x 3):
R = cov(ym);

%   Scale it to suitable form for LSQsolve (3N x 3N):
R = kron(eye(3*N), R);

%   The measurement modelfunction handle (scaled from :
g = @(p) kron(ones(N,1), [eye(3) (3*p*p' - norm(p).^2*eye(3))/ ...
        norm(p)^5]*th_estimate);
    
%   Jacobian function handle
Jacobian = @(p) G_jacobian_new(p, N, th_estimate(4:6));

%   Finally it's time to estimate the location using given function:

[theta, thetas, Js, gammas] = lsqsolve(ym, g, Jacobian, R, theta0);

%   Next we can plot the initial guess and the results

%   Plot of the cost function

    
%   Actual Jacobian function handle:
function G_expanded = G_jacobian_new(p, N, m)
    A = [4*p(1) 3*p(2) 3*p(3); 3*p(2) -2*p(1) 0; 3*p(3) 0 -2*p(1)];
    B = [-2*p(2) 3*p(1) 0; 3*p(1) 4*p(2) 3*p(3); 0 3*p(3) -2*p(2)];
    
    column1 = (A./norm(p)^5 - 5*p(1)*(3*p*p'-norm(p)^2*eye(3))/norm(p)^7)*m;
    column2 = (B./norm(p)^5 - 5*p(2)*(3*p*p'-norm(p)^2*eye(3))/norm(p)^7)*m;
    column3 = zeros(3,1);
    
    G = [column1 column2 column3];
    
    %   Scale the model to be (3N x 2)
    G_expanded = kron(ones(N,1), G);
    
end


%   Testing purposes:
%G_theta = Jacobian(theta0)