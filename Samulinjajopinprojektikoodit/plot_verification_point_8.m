clear
clc

p_verif = [60 0 10; 30 30 10; 45 45 10; 30 45 10; 0 60 10;
           45 30 10; -30 15 10; -30 -15 10]';

y8 = load_data('sLogVal_pos8.txt');         % Measurement data
p8 = [-32 -32 10]';                         % Initial parameter guess

G_symbolic = jacobian_symbolic();

[theta8, thetas8, g, R] = estimator(y8, p8, G_symbolic);

y8 = y8(:);

% Draw plots
J = @(p, y) (y - g(p))'/R*(y - g(p));           % Cost function

pxx = linspace(-35, -20);
pyy = linspace(-35, -20);

Js = zeros(size(pxx));
for i = 1:size(pxx, 2)
   for j = 1:size(pyy, 2)
      Js(i, j) = J([pxx(i); pyy(j); 0], y8);
   end
end
[pxx, pyy] = meshgrid(pxx, pyy);

% [pxx, pyy] = meshgrid(-50:1:70, -50:1:70);      % Grid
% 
% Js = zeros(size(pxx));
% for i = 1:size(pxx, 1)
%     for j = 1:size(pyy, 2)
%         Js(i, j) = J([pxx(i, j); pyy(i, j); 0], y8);
%     end
% end

figure();
contour(pxx, pyy, log(Js), 300);
hold on;
plot(-30, -25, '+k')                                 % Correct position

plot(theta8(1), theta8(2), 'oc');
plot(thetas8(1, :), thetas8(2, :), '-xr');

axis equal;