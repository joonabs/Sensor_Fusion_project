% Static target localization example
%
% 2018-10-11 -- Roland Hostettler <roland.hostettler@aalto.fi>

clear variables;
rng(188);

%% Parameters
% Target parameters: p
pt = [5 3.7]';

% Sensor locations
ps = [
    2, 0, 10, 7;
    6, 0,  2, 8;
];

% Covariances, the ith column corresponds to the ith sensor in ps
sigma2r = [
    1, 0.1, 0.8, 0.5;
    0.3, 0.5, 2, 1.5;
];

%% Model
% Linear model that measures the position of the target; we stack the
% measurements of the sensor on top of each other.
Ns = length(ps);
Gn = eye(2);
G = kron(ones(Ns, 1), Gn);
R = diag(sigma2r(:));

%% Generate data
r = sqrt(R).'*randn(2*Ns, 1);
y = G*pt - ps(:) + r;

%% Estimate
% Least Squares
pthat_ls = G\(y+ps(:));
P_ls = (G'*G)\(G'*R*G)/((G'*G)');

% Weighted least squares
pthat_wls = (G'/R*G)\(G'/R)*(y+ps(:));
P_wls = (G'/R*G)\eye(2);

% Sequential least squares
yy = reshape(y, [2, Ns]);
pthat_sls = zeros(2, Ns);
P_sls = zeros(2, 2, Ns);
pthat_sls(:, 1) = (Gn'/diag(sigma2r(:, 1))*Gn)\(Gn'/diag(sigma2r(:, 1)))*(yy(:, 1)+ps(:, 1));
P_sls(:, :, 1) = (Gn'/diag(sigma2r(:, 1))*Gn)\eye(2);
for i = 2:Ns
    L = P_sls(:, :, i-1)*Gn'/(Gn*P_sls(:, :, i-1)*Gn' + diag(sigma2r(:, i)));
    pthat_sls(:, i) = pthat_sls(:, i-1) + L*(yy(:, i)+ps(:, i) - Gn*pthat_sls(:, i-1));
    P_sls(:, :, i) = P_sls(:, :, i-1) - L*(Gn*P_sls(:, :, i-1)*Gn' + diag(sigma2r(:, i)))*L';
end

% Regularized least squares, case I: mean is close, covariance is large
m0 = [5 5].';
P0 = 2^2*eye(2);

K = P0*G'/(G*P0*G' + R);
pthat_rls = m0 + K*(y+ps(:) - G*m0);
P_rls = P0 - K*(G*P0*G' + R)*K';

% Regularized least squares, case II: mean is far, covariance is small
m0 = [0 0].';
P0 = 1*eye(2);

K = P0*G'/(G*P0*G' + R);
pthat_rls2 = m0 + K*(y+ps(:) - G*m0);
P_rls2 = P0 - K*(G*P0*G' + R)*K';

%% Illustrations
figure(1); clf();

% Sensor & target
plot(ps(1, :), ps(2, :), 'o'); hold on;
plot(pt(1), pt(2), 'x');

% Estimates
colorIndex = get(gca, 'ColorOrderIndex');
plot(pthat_ls(1), pthat_ls(2), 'o');
plot(pthat_wls(1), pthat_wls(2), 'o');
plot(pthat_sls(1, :), pthat_sls(2, :), '--o');
plot(pthat_rls(1), pthat_rls(2), 'o');
plot(pthat_rls2(1), pthat_rls2(2), 'o');
legend('Sensors', 'Target', 'LS', 'WLS', 'SLS', 'RLS (1)', 'RLS (2)');

% Uncertainties
set(gca, 'ColorOrderIndex', colorIndex);
c = cos(linspace(0, 2*pi));
s = sin(linspace(0, 2*pi));
plot(pthat_ls(1) + 2*sqrt(P_ls(1, 1))*c, pthat_ls(2) + 2*sqrt(P_ls(2, 2))*s, '--');
plot(pthat_wls(1) + 2*sqrt(P_wls(1, 1))*c, pthat_wls(2) + 2*sqrt(P_wls(2, 2))*s, '--');
plot(pthat_sls(1, end) + 2*sqrt(P_sls(1, 1, end))*c, pthat_sls(2, end) + 2*sqrt(P_sls(2, 2, end))*s, '--');
plot(pthat_rls(1) + 2*sqrt(P_rls(1, 1))*c, pthat_rls(2) + 2*sqrt(P_rls(2, 2))*s, '--');
plot(pthat_rls2(1) + 2*sqrt(P_rls2(1, 1))*c, pthat_rls2(2) + 2*sqrt(P_rls2(2, 2))*s, '--');
axis equal;