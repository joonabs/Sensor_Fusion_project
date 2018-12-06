function [theta, thetas, Js, gammas] = lsqsolve(y, g, G, R, theta0, Imax, method)
% Nonlinear least squares solver
%
% USAGE
%   theta = LSQSOLVE(y, g, G, R, theta0)
%   [theta, thetas, Js, gammas] = LSQSOLVE(y, g, G, R, theta0, Imax, method)
%
% DESCRIPTION
%   Nonlinear least squares solver implementing the gradient dscent,
%   Gauss-Newton, and Levenberg-Marquardt algorithms. The default method is 
%   gradient descent.
%
% PARAMETERS
%   y       Measurement data
%   g       Measurement model, function handsize(le g(theta)
%   G       Jacobian matrix of the measurement model, function handle 
%           G(theta)
%   R       Measurement noise covariance matrix
%   theta0  Initial parameter guess
%   Imax    Maximum number of iterations (optional, default: 10)
%   method  Method to use, one of 'gradient', 'gauss-newton', or
%           'levenberg-marquardt' (optional, default: 'gradient')
%
% RETURNS
%   theta   Parameter estimate
%   thetas  Parameters for each iteration
%   Js      Cost for each iteration
%   gammas  Step lengths (for gradient descent and Gauss-Newton) or damping
%           factor (for Levenberg-Marquardt)
%
% AUTHOR
%   2018-09-28 -- Roland Hostettler <roland.hostettler@aalto.fi>

    %% Defaults
    narginchk(5, 7);
    if nargin < 6 || isempty(Imax)
        Imax = 10;
    end
    if nargin < 7 || isempty(method)
        method = 'gradient';
    end
    
    %% Parameters
    % Stopping tolerance
    epsilon = 1e-3;
        
    % Parameters for Levenberg-Marquardt (used for adapting the damping)
    nu = 2;
    tau = 1;
    
    % Step length (gradient descent and Gauss-Newton) or damping constant
    % (Levenberg-Marquardt)
    gamma = tau*max(diag(G(theta0)'/R*G(theta0)));
    
    %% Preallocate
    Ntheta = size(theta0, 1);
    thetas = zeros(Ntheta, Imax+1);
    Js = zeros(1, Imax+1);
    gammas = zeros(1, Imax);

    %% Initialize
    % Cost function
    LR = chol(R).';
    Jwls = @(theta) sum((LR\(y-g(theta))).^2);

    % Initial parameters and cost
    theta = theta0;
    J = Jwls(theta0);
    Js(1) = J;
    thetas(:, 1) = theta0;

    
    %% Iterations
    done = false;
    i = 1;
    while ~done
        %% Calculate the descent direction
        if strcmpi(method, 'gradient')
            dtheta = G(theta)'/R*(y - g(theta));
        elseif strcmpi(method, 'gauss-newton')
            GG = G(theta);
            dtheta = (GG'/R*GG)\(GG'/R)*(y - g(theta));
        elseif strcmpi(method, 'levenberg-marquardt')
            GG = G(theta);
            dtheta = @(lambda) (GG'/R*GG + lambda*eye(Ntheta))\(GG'/R)*(y - g(theta));
        else
            error('Unknown method.');
        end

        %% Parameter update
        if strcmp(method, 'levenberg-marquardt')
            %% Levenberg-Marquardt
            done2 = false;
            while ~done2
                dthetap = dtheta(gamma);
                thetap = theta + dthetap;
                Jp = Jwls(thetap);
                
                rho = (J-Jp)/(2*dthetap'*GG'/R*(y-g(theta)) - dthetap'*GG'/R*GG*dthetap);
                if rho >= 0
                    % Decrease in cost, decrease damping
                    nu = 2;
                    gamma = tau*gamma*max(1/3, 1-(2*rho-1)^3);
                    done2 = true;
                else
                    % Increase in cost, increase damping
                    gamma = nu*gamma;
                    nu = 2*nu;
                    warning('Damping increased');
                end
            end 
        else
            %% Gradient descent and Gauss-Newton
            gamma = 1;
            
            done2 = false;
            while ~done2
                gamma = gamma/2;
                    thetap = theta + gamma*dtheta;

                Jp = Jwls(thetap);
                done2 = (Jp <= J);
            end
        end
        
        %% Convergence check
        done = ( ...
            i >= Imax ...                                           % Maximum number of iterations
            || Jp-J == 0 ...                                        % No decrease in cost
            || abs(Jp-J)/J < epsilon ...                            % Relative decrease in cost
            || abs(norm(thetap-theta))/norm(theta) < epsilon ...    % Relative decrease in parameter magnitude
        );

        
        %% Parameter update
        theta = thetap;
        J = Jp;
        i = i+1;
        thetas(:, i) = theta;
        Js(i) = J;
        gammas(i) = gamma;        
    end
    
    %% Strip extra entries
    thetas = thetas(:, 1:i);
    Js = Js(1:i);
    gammas = gammas(1:i);
end
