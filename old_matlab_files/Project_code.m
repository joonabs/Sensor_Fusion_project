file1 = "sensorLog_20181025T135642.txt";
file2 = "sensorLog_20181025T135730.txt";
[ym1, tm1, ya1, ta1] = load_data(file1);
[ym2, tm2, ya2, ta2] = load_data(file2);
p = [30 0 10 ; 30 30 10]';  % aleksin paperin mittauspisteet!!
H = (3*p*p'- norm(p).^2*eye(3))/(norm(p).^5);  
%%
G = [eye(3) H];
G1 = [repmat(G,size(ym1)/3),1];
G2 = [repmat(G,size(ym2)/3),1];
G = [G1(:) ; G2(:)];
Theta_hat_least_squares = (G'*G)\G'*ym;



%%
theta0 = [1 2 3]'; % alkuarvaus lähelle oikeaa pistettä! 
y = ym;
g = ; % jatka tästä
G = %@(pt) (1./g(pt)*ones(1, 2)).*(pt*ones(1, Ns)-ps).'; % Jacobian
R = cov(y')'; % Kumminpäin tämä eli tuleeko tupla transponointi
Imax = 1000;
method = 'levenberg-marquardt';
[theta, thetas, Js, gammas] = lsqsolve(y, g, G, R, theta0, Imax, method);