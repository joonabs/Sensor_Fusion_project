clear all
clc

%% Calibrate
y = load_data('Verification_data_5.txt');       % Measurement data
p0 = [20 -39]';                                 % Initial parameter guess
target = [20 -30]';                             % Target location
z = 11;
G_symbolic = jacobian_symbolic();
N = size(y,2);                                  % Number of measurements

th_estimate = calibration();

% Measurement model, function handle g(theta)
g = @(p) kron(ones(N,1), [eye(3) (3*[p; z]*[p; z]'-norm([p; z]).^2*eye(3))/norm([p; z]).^5]*th_estimate);

R = kron(eye(N), cov(y'));                      % Measurement noise covariance matrix    
Imax = 20;                                      % Maximum number of iterations (optional, default: 10)    
y = y(:);

%% Draw cost function
J = @(p, y) (y - g(p))'/R*(y - g(p));           % Cost function

[pxx, pyy] = meshgrid(-10:1:50, -50:1:10);      % Grid

Js = zeros(size(pxx));
for i = 1:size(pxx, 1)
    for j = 1:size(pyy, 2)
        Js(i, j) = J([pxx(i, j); pyy(i, j)], y);
    end
end

figure();
hold on
axis equal
contour(pxx, pyy, log(Js), 50, 'DisplayName', 'Cost function');
plot(0, 0, '+r', 'DisplayName', 'Sensor location: (0, 0)');
plot(target(1), target(2), '+k', 'DisplayName', ['True location: (', num2str(target(1)), ', ', num2str(target(2)), ')'], 'LineWidth', 2);


%% Calculate and draw optimization
methods = {'gradient', 'gauss-newton', 'levenberg-marquardt'};         
plot_line_styles = {'-xr', '-xk', '-xb'};
names = {'Gradient descent', 'Gauss-Newton', 'Levenberg-Marquardt'};
final_errors = [];

for i = 1:3
    % Jacobian matrix of the measurement model, function handle                 
    G_jacobian = @(p) jacobian_substitute(G_symbolic, p, th_estimate(4:6), N);

    [theta, thetas, Js, gammas] = lsqsolve(y, g, G_jacobian, R, p0, Imax, methods{i});
    final_errors = [final_errors norm(theta - target)];

    plot(theta(1), theta(2), 'oc', 'HandleVisibility','off');
    plot(thetas(1, :), thetas(2, :), plot_line_styles{i}, 'DisplayName', names{i} + ", " + length(gammas) + " iterations");
end
final_error = round(mean(final_errors),1);
title(['Mean error from true location: ', num2str(final_error), ' cm'], 'FontSize', 14);
legend('Location', 'southoutside', 'FontSize', 14)