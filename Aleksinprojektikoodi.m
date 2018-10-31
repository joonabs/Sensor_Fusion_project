% Project work - static part
% Aleksi Mäkinen and Joonas Isometsä
clear
clc
%   The initial parameters required for the code to work:
%       
%   Measurements to calibrate the model from measurements-file, at least 2 needed!
file11 = 'sensorLog_20181025T135642.txt';
file12 = 'sensorLog_20181025T135730.txt';
%file11 = 'sensorLog_20181025T135926.txt';
%file12 = 'sensorLog_20181025T135823.txt';

%   The corresponding known actual positions:
%   [x1, y1, z1; x2, y2, z2; ... etc]
p = [30 0 10; 30 30 10]';
%p = [30 -30 10; 0 45 10]';


%   Measurements to validate the model from measurements-file, at least 2 needed!
file21 = 'sensorLog_20181025T135823.txt';
file22 = 'sensorLog_20181025T135926.txt';
%   Initial position guess for validation:
%   [x0, y0, z0]'
theta0 = [10,10,10];
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
[ym1, tm, ya, ta] = load_data(file11);
[ym2, tm, ya, ta] = load_data(file12);

%   Add the measurements of the 2 positions to the same 2*L x 1 matrix
y = [ym1(:); ym2(:)];
%   The known positions p are used here.
H = (3*p*p' - norm(p).^2*eye(3)) / norm(p).^5;
G = [eye(3) H];
%   Generate appropriately sized G to match the y measurements
G = repmat(G, size(y,1)/3, 1);

% Now the estimate can be calculated:
th_estimate = (G'*G)\G'*y


%       Then the validation
%
%       The following values and functions are needed:
%
%   y       Measurement data
%   g       Measurement model, function handle g(theta)
%   G       Jacobian matrix of the measurement model, function handle 
%           G(theta)
%   R       Measurement noise covariance matrix
%   theta0  Already set
%   Imax    - || -
%   method  - || -
%       
%   We start by loading the validation data:
[ym1, tm, ya, ta] = load_data(file21);
[ym2, tm, ya, ta] = load_data(file22);
%   Then we call the measurement modelfunction:
function [g_theta] = g(theta)

end

%Testing purposes:
%[G_theta] = Jacobian(theta0)

% Jacobian function handle:
function[G_theta] = Jacobian(theta0)
    syms p_x p_y p_z
    p = [p_x; p_y; p_z];
    M = (3*p*p' - sqrt(p_x^2 + p_y^2 + p_z^2)^2)...
        /sqrt(p_x^2 + p_y^2 + p_z^2)^5;
    J_x = jacobian(M(1,:), [p_x, p_y, p_z]);
    J_y = jacobian(M(2,:), [p_x, p_y, p_z]);
    J_z = jacobian(M(3,:), [p_x, p_y, p_z]);
    G_theta(p_x, p_y, p_z) = [J_x; J_y; J_z];
    G_theta = G_theta(theta0(1), theta0(2), theta0(3))
end
%p_trans = [p_x p_y p_z];

% Now the measurement noise covariance matrix R:
%R = cov(