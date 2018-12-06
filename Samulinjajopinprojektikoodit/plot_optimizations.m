clear
clc

p_verif = [60 0 10; 30 30 10; 45 45 10; 30 45 10; 0 60 10; 45 30 10; -30 15 10; -30 -15 10]';

y7 = load_data('sLogVal_pos7.txt');         % Measurement data
p7_1 = [70 20 10]';                         % Initial parameter guess
p7_2 = [0 60 10]';
p7_3 = [-40 -40 10]';
p7_4 = [60 -30 10]';

G_symbolic = jacobian_symbolic();

[theta7_1, thetas7_1, g, R] = estimator(y7, p7_1, G_symbolic);
[theta7_2, thetas7_2] = estimator(y7, p7_2, G_symbolic);
[theta7_3, thetas7_3] = estimator(y7, p7_3, G_symbolic);
[theta7_4, thetas7_4] = estimator(y7, p7_4, G_symbolic);

theta7 = {theta7_1, theta7_2, theta7_3, theta7_4};
thetas7 = {thetas7_1, thetas7_2, thetas7_3, thetas7_4};
y7 = y7(:);

% Draw plots
J = @(p, y) (y - g(p))'/R*(y - g(p));           % Cost function
[pxx, pyy] = meshgrid(-50:1:70, -50:1:70);      % Grid

Js = zeros(size(pxx));
for i = 1:size(pxx, 1)
    for j = 1:size(pyy, 2)
        Js(i, j) = J([pxx(i, j); pyy(i, j); 0], y7);
    end
end

figure();
contour(pxx, pyy, log(Js), 500);
%title('Four different initial guesses - Measurement point [-30 15]')
hold on;
plot(-30, 15, '+k')                                 % Correct position

for i = 1:4
    plot(theta7{i}(1), theta7{i}(2), 'oc');
    plot(thetas7{i}(1, :), thetas7{i}(2, :), '-xr');
end

axis equal;


%%
clear
clc

p_verif = [60 0 10; 30 30 10; 45 45 10; 30 45 10;
    0 60 10; 45 30 10; -30 15 10; -30 -15 10]';

p = [50 -10 10; 40 40 10; 55 55 10; 40 50 10]';  % Initial parameter guesses

G_symbolic = jacobian_symbolic();

figure();
[pxx, pyy] = meshgrid(-50:1:70, -50:1:70);      % Grid
Js = zeros(size(pxx));

i = 1;
for i = 1:4
    subplot(2,2,i);
    filename = sprintf('sLogVal_pos%i.txt',i);
    y = load_data(filename);                  % Measurement data

    [theta, thetas, g, R] = estimator(y, p(:,i), G_symbolic);
    J = @(p, y) (y - g(p))'/R*(y - g(p));           % Cost function

    y = y(:);

    for j = 1:size(pxx, 1)
        for k = 1:size(pyy, 2)
            Js(j, k) = J([pxx(j, k); pyy(j, k); 0], y);
        end
    end
    contour(pxx, pyy, log(Js), 500);
    hold on


    plot(p_verif(1,i), p_verif(2,i), '+c')                                 % Correct position
    plot(theta(1), theta(2), 'ok');
    plot(thetas(1, :), thetas(2, :), '-xr');
    axis equal
end
