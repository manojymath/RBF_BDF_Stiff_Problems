
% % % % "RBF based backward differentiation methods for stiff differential equations"
% % % % Authors: A. Sreedhar, Manoj Kumar Yadav, Chirala Satyanarayana


% % % % % Define the ODEs
% % % % % Example 1: dy/dt = lambda*(u-cos(t)); w.r.to u(0)=1;
% % % % % Example 2: dy/dt = lambda*(u-cos(t))-sin(t); w.r.to u(0)=1.5;
% % % % % To solve by using the BDF2, GA-BDF2 and MQ-BDF2 methods;
clc;
clear all;

% % % % % Define the Stiff ODE for IVP (BDF-Paper-example 1)
lbda = -10^(4); % % % % Parameter Value
f = @(u, t) lbda*(u-cos(t));  % % % % % Given function
uexact = @(t) (lbda^2/(1+lbda^2))*cos(t)-(lbda/(1+lbda^2))*sin(t)+...
              (1/(1+lbda^2))*exp(lbda*t); % % % %  Exact Solution
a = 0; b = 7;      % % % % % Time interval
u0 = 1; u(1) = u0; % % % % % Initial condition

% % % % % Define the Stiff ODE for IVP (BDF-Paper-example 2)
% lambda = -10^(6);
% f = @(u, t) lambda*(u-cos(t))-sin(t);    % % % % Given function
% uexact = @(t) exp(lambda*t)*(0.5)+cos(t);% % % % Exact Solution
% a = 0; b = 3;        % % % % Time interval
% u0 = 1.5; u(1) = u0; % % % % Initial condition

% % % % % Array of M values
M_values = [11, 21, 41, 81, 161, 321, 641];
errors = zeros(size(M_values)); % % % % Intialization for Errors

for k = 1:length(M_values)
    M = M_values(k);
    t = linspace(a, b, M);
    dt = (b-a)/(M-1);

    % % % % Estimating the u(2) solution
    u(2) = uexact(t(2));

    % % % Iterate using the provided weights and respective BDF2 method;
    % % % % %  Solving the ODE using BDF2 method
    for i = 3:M
        ucurrent = u(i-1); uprevious = u(i-2);% % Initial guesses for u_n
        % % % % % Define BDF2 weights
        alphac1 = 3/(2*dt); alphac2 = -2/dt; alphac3 = 1/(2*dt);
        % % % % % Define the function to find the root
        F2 = @(un) alphac1*un+alphac2*ucurrent+alphac3*uprevious-f(un,t(i));
        options = optimset('Display', 'off'); % % Suppress output of fsolve
        % % % % % Use fsolve to find the next value
        u(i) = fsolve(F2, ucurrent, options);
    end

    % % % % % Call the derivative values by Fornberg algorithm
    [u1,u2,u3] = FD_deri_values(t, u);

    % % % % % % For the optimised epsilon value and sol'n Intialisation
    ep=zeros(1,M); uRBDF=zeros(1,M);

    % % % To intialize the first solution
    uRBDF(1) = u(1); ucurrent = uRBDF(1); % % % % % Intial condition

    % % % % Estimating the u(2) solution
    uRBDF(2) = uexact(t(2));

    % % % % % Now apply RBF-BDF2 for subsequent steps
    for i=3:M
        ucurrent = uRBDF(i-1); uprevious = uRBDF(i-2);

        % % % % % Calculate epsilon(i)
        % ep(i) = -(u3(i))/(6*u1(i)); % % % % opt. epsilon with GA-Weights
        ep(i) = -(u3(i))/(3*u1(i)); % % % % opt. epsilon with MQ-Weights

        % % % % % Define GA-BDF2 weights
        % alpha1 = (3/(2*dt))-(11*ep(i)*dt)/6+(37*ep(i)^2*dt^3)/36+(445*ep(i)^3*dt^5)/216-(2129*ep(i)^4*dt^7)/1296;
        % alpha2 = (-2/dt)+(5*ep(i)*dt)/3-(13*ep(i)^2*dt^3)/18-(445*ep(i)^3*dt^5)/108-(3631*ep(i)^4*dt^7)/648;
        % alpha3 = (1/(2*dt))+(ep(i)*dt)/6-(11*ep(i)^2*dt^3)/36+(445*ep(i)^3*dt^5)/216+(9391*ep(i)^4*dt^7)/1296;

        % % % % Define MQ-BDF2 weights
        alpha1 = (3/(2*dt))-(7*ep(i)*dt)/4+(79*ep(i)^2*dt^3)/16+(755*ep(i)^3*dt^5)/64;
        alpha2 = (-2/dt)+(5*ep(i)*dt)/2-(55*ep(i)^2*dt^3)/8-(1075*ep(i)^3*dt^5)/32;
        alpha3 = (1/(2*dt))-(3*ep(i)*dt)/4+(31*ep(i)^2*dt^3)/16+(1395*ep(i)^3*dt^5)/64;

        % % % % % Define the function to find the root
        F4 = @(un) alpha1*un+alpha2*ucurrent+alpha3*uprevious-f(un,t(i));
        options = optimset('Display', 'off'); % Suppress output of fsolve
        % % % % % Use fsolve to find the next value
        uRBDF(i) = fsolve(F4, ucurrent, options);
    end
    % % % % % % Calculating the Error
    errors(k) = abs(uexact(t(M))-uRBDF(M));
end

% % % % % % Calculate the order of convergences
orders = size(errors); % Initialize with zeros
orders(1) = 0;
for j = 2:length(M_values)
    orders(j) = log(errors(j-1)/errors(j))/log(2);
end

% % % % % Display results in scientific notation
errors_scientific = sprintfc('%0.4e', errors);
table(M_values', errors_scientific', orders', 'VariableNames', {'M', 'Error', 'Order'})
