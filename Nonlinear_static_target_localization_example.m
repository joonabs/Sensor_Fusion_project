% Nonlinear static target localization example
%
% 2018-10-12 -- Roland Hostettler <roland.hostettler@aalto.fi>

clear variables;
%addpath lib;
rng(188);

%% Parameters
% True target location
pt = [5 3.7]';

% Initial guess of the target location
% pt0 = [1, 2].';     % First example, close to true position
pt0 = [-5, 0].';    % Second example, farther away

% Sensor locations, one column per sensor [px, py]^T
ps = [
    2, 0, 10, 7;
    6, 0,  2, 8;
];

% Covariances, one column per sensor
sigma2r = [
    1, 0.1, 0.8, 0.5;
];

% Maximum number of iterations for the numerical solvers
Imax = 20;

%% Model
Ns = length(ps);
g = @(pt) sqrt(sum((pt*ones(1, Ns) - ps).^2, 1)).';     % Measurement model (range)
G = @(pt) (1./g(pt)*ones(1, 2)).*(pt*ones(1, Ns)-ps).'; % Jacobian
R = diag(sigma2r(:));

%% Generate data
r = sqrt(R).'*randn(Ns, 1);
y = g(pt) + r;

%% Estimation
% Gradient descent
[pthat_gd, pts_gd, Js_gd] = lsqsolve(y, g, G, R, pt0, Imax, 'gradient');

% Gauss-Newton
[pthat_gn, pts_gn, Js_gn] = lsqsolve(y, g, G, R, pt0, Imax, 'gauss-newton');

% Levenberg-Marquardt
[pthat_lm, pts_lm, Js_lm] = lsqsolve(y, g, G, R, pt0, Imax, 'levenberg-marquardt');

%% Calculate contours
J = @(pt) (y - g(pt))'/R*(y-g(pt));                     % Cost function
% [pxx, pyy] = meshgrid(0:0.25:10, 0:0.25:10);            % Grid for first example
% [pxx, pyy] = meshgrid(4:0.1:6, 3:0.1:5);                % Zoomed grid
[pxx, pyy] = meshgrid(-7.5:0.25:7.5, -2.5:0.25:12.5);   % Grid for second example
Js = zeros(size(pxx));
for iu = 1:size(pxx, 1)
    for iv = 1:size(pyy, 2)
        Js(iu, iv) = J([pxx(iu, iv); pyy(iu, iv)]);
    end
end

%% Illustrations
figure(1); clf();
contour(pxx, pyy, Js, 20); hold on;
plot(ps(1, :), ps(2, :), 'ok');
plot(pt(1), pt(2), 'xr');

plot(pthat_gd(1), pthat_gd(2), 'ob');
plot(pthat_gn(1), pthat_gn(2), 'or');
plot(pthat_lm(1), pthat_lm(2), 'ok');

plot(pts_gd(1, :), pts_gd(2, :), '-xb');
plot(pts_gn(1, :), pts_gn(2, :), '-xr');
plot(pts_lm(1, :), pts_lm(2, :), '-xk');

axis([0, 10, 0, 10]);
axis equal;