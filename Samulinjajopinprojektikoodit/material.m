clear
clc
syms px py pz
syms m1 m2 m3
p = [px; py; pz];
p_trans = [px py pz];

m = [m1 m2 m3]';
M = ((3*p*p_trans - sqrt(px^2 + py^2 + pz^2)^2)/sqrt(px^2 + py^2 + pz^2)^5)
M*m
size(M)
%A = jacobian(M, [px, py, pz])

%%

A = jacobian(M(1,:), [px, py, pz]);
B = jacobian(M(2,:), [px, py, pz]);
C = jacobian(M(3,:), [px, py, pz]);

G = [A; B; C]

G = subs(G, px, 1)
G = subs(G, py, 2)
G = subs(G, pz, 3)



%%
syms r l f 
x = r*cos(l)*cos(f); y = r*cos(l)*sin(f); z = r*sin(l);
J = jacobian([x; y; z], [r l f])

%%
% Jacobian
syms px py pz
syms mx my mz
p = [px; py; pz];
p_trans = [px py pz];
m = [mx; my; mz];

%M = ((3*p*p_trans - sqrt(px^2 + py^2 + pz^2)^2)/sqrt(px^2 + py^2 + pz^2)^5)*th_estimate(4:6);
M = ((3*p*p_trans - sqrt(px^2 + py^2 + pz^2)^2)/sqrt(px^2 + py^2 + pz^2)^5)*m;

A = simplify(jacobian(M, px))
%%
A = jacobian(M, [px, py]);
A = [A 0]
simplify(A)

G_jacobian = @(p) kron(ones(N,1), subs(A, {px, py, pz}, {p(1), p(2), p(3)}));

%%
%     A = [(4*m(1)*p(1))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (m(2)*(2*p(1) - 3*p(2)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (m(3)*(2*p(1) - 3*p(3)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (5*m(2)*p(1)*(p(1)^2 - 3*p(1)*p(2) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(3)*p(1)*(p(1)^2 - 3*p(1)*p(3) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(1)*p(1)*(- 2*p(1)^2 + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2),        (m(2)*(3*p(1) - 2*p(2)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (2*m(3)*p(2))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (2*m(1)*p(2))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (5*m(2)*p(2)*(p(1)^2 - 3*p(1)*p(2) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(3)*p(2)*(p(1)^2 - 3*p(1)*p(3) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(1)*p(2)*(- 2*p(1)^2 + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2), 0];
%     B = [(5*m(1)*p(1)*(p(1)^2 - 3*p(1)*p(2) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) - (2*m(3)*p(1))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (m(1)*(2*p(1) - 3*p(2)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (2*m(2)*p(1))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (5*m(3)*p(1)*(p(1)^2 + p(2)^2 - 3*p(2)*p(3) + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(2)*p(1)*(p(1)^2 - 2*p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2), (4*m(2)*p(2))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (m(1)*(3*p(1) - 2*p(2)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (m(3)*(2*p(2) - 3*p(3)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (5*m(1)*p(2)*(p(1)^2 - 3*p(1)*p(2) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(3)*p(2)*(p(1)^2 + p(2)^2 - 3*p(2)*p(3) + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(2)*p(2)*(p(1)^2 - 2*p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2), 0];
%     C = [(5*m(1)*p(1)*(p(1)^2 - 3*p(1)*p(3) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) - (2*m(3)*p(1))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (m(1)*(2*p(1) - 3*p(3)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (2*m(2)*p(1))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (5*m(2)*p(1)*(p(1)^2 + p(2)^2 - 3*p(2)*p(3) + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(3)*p(1)*(p(1)^2 + p(2)^2 - 2*p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2),          (5*m(1)*p(2)*(p(1)^2 - 3*p(1)*p(3) + p(2)^2 + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) - (2*m(3)*p(2))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (m(2)*(2*p(2) - 3*p(3)))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) - (2*m(1)*p(2))/(p(1)^2 + p(2)^2 + p(3)^2)^(5/2) + (5*m(2)*p(2)*(p(1)^2 + p(2)^2 - 3*p(2)*p(3) + p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2) + (5*m(3)*p(2)*(p(1)^2 + p(2)^2 - 2*p(3)^2))/(p(1)^2 + p(2)^2 + p(3)^2)^(7/2), 0];
% 
%     G = [A; B; C];
%     
%     G_expanded = kron(ones(N,1), G);

%%
% th1grid = linspace(-60,60);
% th2grid = linspace(-60,60);
% costFunc = zeros(length(th1grid), length(th2grid));
% 
% for i = 1:length(th1grid)
%     for j = 1:length(th2grid)
%         thTemp = [th1grid(i); th2grid(j); 0];
%         costFunc(i,j) = J(thTemp);
%     end
% end

%contour(th1grid, th2grid, costFunc, 100);  % from exercise 3.1
