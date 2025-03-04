% % % % "RBF based backward differentiation methods for stiff differential equations"
% % % % Authors: A. Sreedhar, Manoj Kumar Yadav, Chirala Satyanarayana


% % % To solve the system of stiff equations using BDF methods
% % % Define the Sytem of ODEs of the form du/dt=f(t,u(t),v(t)); and dv/dt=g(t,u(t),v(t));
% % % w.r.to the conditions u(t_0)=u0; and v(t_0)=v0; 
clc
clear all

% % % % Define the system ODEs (Example-4)
% % % % Define the intialization
a=0; b=100; u0 = 2; v0 = 2;  % % % % Intial conditions
% % % % % % Define the exact solutions
uexact = @(t) exp(-1000*t).*(sin(1000*t)+cos(1000*t))+sin(t)+exp(-t); % % Exact Solution for u
vexact = @(t) exp(-1000*t).*(cos(1000*t)-sin(1000*t))+cos(t); % % Exact Solution for v

% % % % % Array of M values for Number of time steps (nodes)
M_values = [321, 641, 1281, 2561, 5121];

% % % Intialization of errors Matrices w.r.to the solutions 'u'and 'v'
errorsu1 = zeros(size(M_values)); errorsv1 = zeros(size(M_values));
errorsu2 = zeros(size(M_values)); errorsv2 = zeros(size(M_values));
errorsu3 = zeros(size(M_values)); errorsv3 = zeros(size(M_values));

for k = 1:length(M_values)
    M = M_values(k);
    t = linspace(a, b, M); % % % % Define the time grid
    dt = (b-a)/(M-1); % % % % Step size, calculated from no.of nodes

    % % % % Preallocate the solution arrays
    P1 = zeros(2, M); % % % % (for BDF)
    P2 = zeros(2, M); % % % % (For GA-BDF)
    P3 = zeros(2, M); % % % % (For MQ-BDF)

    % % % % Setting Initial conditions
    P0 = [u0; v0];
    P1(:,1) = P0; % % % % Set initial conditions at t = 0 for BDF
    P2(:,1) = P0; % % % % Set initial conditions at t = 0 for GA-BDF
    P3(:,1) = P0; % % % % Set initial conditions at t = 0 for MQ-BDF

    % % % % Find the required Intial solution using BDF1 method
    % % % % Define BDF1 weights
    alphaB11 = 1/dt; alphaB10 = -1/dt;
     betaB11 = 1/dt;  betaB10 = -1/dt;

    % % % % To define the required Matrices
    A11 = [alphaB11 0; 0 betaB11];  B11 = [alphaB10 0; 0 betaB10];
    D11 = [-1000 1000; -1000 -1000];
    F11 = [999*exp(-t(2))+1000*sin(t(2))-999*cos(t(2)); 1000*exp(-t(2))+999*sin(t(2))+1000*cos(t(2))];
    G11 = A11-D11; X11 = -(G11\B11);

    % % % % Scheme to find the solution
    P1(:, 2) = X11*P1(:,1)+(G11\F11);

    % % % % % % Apply Classical BDF2 for remaining steps
    for i=3:M

        % % % % % % Define BDF2 weights
        alphac22 = 3/(2*dt); alphac21 = -2/dt;  alphac20 = 1/(2*dt);
        betac22  = 3/(2*dt); betac21  = -2/dt;  betac20  = 1/(2*dt);

        % % % % To find the required Matrices
        A12 = [alphac22 0; 0 betac22];  B12 = [alphac21 0; 0 betac21]; C12 = [alphac20 0; 0 betac20];
        D12 = [-1000 1000; -1000 -1000];  % % % For Sci. Compu.Paper
        F12 = [999*exp(-t(i))+1000*sin(t(i))-999*cos(t(i)); 1000*exp(-t(i))+999*sin(t(i))+1000*cos(t(i))];
        G12 = A12-D12; X12 = -(G12\B12); Y12 = -(G12\C12);
        % % % % Scheme to find the solution
        P1(:, i) = X12*P1(:,i-1)+Y12*P1(:,i-2)+(G12\F12); % % Scheme for Sci.comp. Paper
    end

    % % % % % % Extract results
    u1 = P1(1,:); v1 = P1(2,:);
    
    % % % % % Call the derivative values by Fornberg algorithm
    [D1u, D2u, D3u] = FD_deri_valuesu1(t, u1);
    [D1v, D2v, D3v] = FD_deri_valuesv1(t, v1);

    % % % For choosing the optimised Epsilon Value Intialisation
    epg1=zeros(1,M); epg2=zeros(1,M); % % % % (For GA-BDF)
    epm1=zeros(1,M); epm2=zeros(1,M); % % % % (For MQ-BDF)

    P2(:,2) = P1(:,2); % % % % Using the classical BDF1 soln
    P3(:,2) = P1(:,2); % % % % Using the classical BDF1 soln

    % % % % % Now apply RBF-BDF2 for subsequent steps
    for i=3:M

        % % % % % Calculate epsilon(i) for GA-BDF2
        epg1(i) = -(D3u(i))/(6*D1u(i));% % % % opt. epsilon with GA-Weights
        epg2(i) = -(D3v(i))/(6*D1v(i));% % % % opt. epsilon with GA-Weights

        % % % % % % Define GA-BDF2 u_weights
        alphag22 = (3/(2*dt))-(11*epg1(i)*dt)/6+(37*epg1(i)^2*dt^3)/36+(445*epg1(i)^3*dt^5)/216-(2129*epg1(i)^4*dt^7)/1296;
        alphag21 = (-2/dt)+(5*epg1(i)*dt)/3-(13*epg1(i)^2*dt^3)/18-(445*epg1(i)^3*dt^5)/108-(3631*epg1(i)^4*dt^7)/648;
        alphag20 = (1/(2*dt))+(epg1(i)*dt)/6-(11*epg1(i)^2*dt^3)/36+(445*epg1(i)^3*dt^5)/216+(9391*epg1(i)^4*dt^7)/1296;

        % % % % Define GA-BDF2 v_weights
        betag22 = (3/(2*dt))-(11*epg2(i)*dt)/6+(37*epg2(i)^2*dt^3)/36+(445*epg2(i)^3*dt^5)/216-(2129*epg2(i)^4*dt^7)/1296;
        betag21 = (-2/dt)+(5*epg2(i)*dt)/3-(13*epg2(i)^2*dt^3)/18-(445*epg2(i)^3*dt^5)/108-(3631*epg2(i)^4*dt^7)/648;
        betag20 = (1/(2*dt))+(epg2(i)*dt)/6-(11*epg2(i)^2*dt^3)/36+(445*epg2(i)^3*dt^5)/216+(9391*epg2(i)^4*dt^7)/1296;


        % % % % % Calculate epsilon(i) for MQ-BDF2
        epm1(i) = -(D3u(i))/(3*D1u(i));% % % % opt. epsilon with MQ-Weights
        epm2(i) = -(D3v(i))/(3*D1v(i));% % % % opt. epsilon with MQ-Weights

        % % % % % Define MQ-BDF2 u_weights
        alpham22 = (3/(2*dt))-(7*epm1(i)*dt)/4+(79*epm1(i)^2*dt^3)/16+(755*epm1(i)^3*dt^5)/64;
        alpham21 = (-2/dt)+(5*epm1(i)*dt)/2-(55*epm1(i)^2*dt^3)/8-(1075*epm1(i)^3*dt^5)/32;
        alpham20 = (1/(2*dt))-(3*epm1(i)*dt)/4+(31*epm1(i)^2*dt^3)/16+(1395*epm1(i)^3*dt^5)/64;

        % % % % Define MQ-BDF2 v_weights
        betam22 = (3/(2*dt))-(7*epm2(i)*dt)/4+(79*epm2(i)^2*dt^3)/16+(755*epm2(i)^3*dt^5)/64;
        betam21 = (-2/dt)+(5*epm2(i)*dt)/2-(55*epm2(i)^2*dt^3)/8-(1075*epm2(i)^3*dt^5)/32;
        betam20 = (1/(2*dt))-(3*epm2(i)*dt)/4+(31*epm2(i)^2*dt^3)/16+(1395*epm2(i)^3*dt^5)/64;

        % % % % To find the required Matrices
        A2 = [alphag22 0; 0 betag22];  B2 = [alphag21 0; 0 betag21]; C2 = [alphag20 0; 0 betag20];
        A3 = [alpham22 0; 0 betam22];  B3 = [alpham21 0; 0 betam21]; C3 = [alpham20 0; 0 betam20];
        D1 = [-1000 1000; -1000 -1000];  % % % For Sci. Compu.Paper
        F = [999*exp(-t(i))+1000*sin(t(i))-999*cos(t(i)); 1000*exp(-t(i))+999*sin(t(i))+1000*cos(t(i))];
        G2 = A2-D1; X2 = -(G2\B2); Y2 = -(G2\C2);
        G3 = A3-D1; X3 = -(G3\B3); Y3 = -(G3\C3);
        % % % % Scheme to find the solution
        P2(:, i) = X2*P2(:,i-1)+Y2*P2(:,i-2)+(G2\F); % % GA-BDF2 scheme solution
        P3(:, i) = X3*P3(:,i-1)+Y3*P3(:,i-2)+(G3\F); % % MQ-BDF2 scheme solution
    end
    % % % % % % Extract results
    u2 = P2(1,:); v2 = P2(2,:); % % % For GA-BDF2
    u3 = P3(1,:); v3 = P3(2,:); % % % For MQ-BDF2

    % % % % % Calculate final time errors
    errorsu1(k) = abs(uexact(t(M))-u1(M));
    errorsv1(k) = abs(vexact(t(M))-v1(M));
    errorsu2(k) = abs(uexact(t(M))-u2(M));
    errorsv2(k) = abs(vexact(t(M))-v2(M));
    errorsu3(k) = abs(uexact(t(M))-u3(M));
    errorsv3(k) = abs(vexact(t(M))-v3(M));
end

% % % Calculate the order of convergence
ordersu1 = log(errorsu1(1:end-1)./errorsu1(2:end))/log(2);
ordersv1 = log(errorsv1(1:end-1)./errorsv1(2:end))/log(2);
ordersu2 = log(errorsu2(1:end-1)./errorsu2(2:end))/log(2);
ordersv2 = log(errorsv2(1:end-1)./errorsv2(2:end))/log(2);
ordersu3 = log(errorsu3(1:end-1)./errorsu3(2:end))/log(2);
ordersv3 = log(errorsv3(1:end-1)./errorsv3(2:end))/log(2);

% % % % Display results in scientific notations
errors_scientificu1 = sprintfc('%0.4e', errorsu1);
errors_scientificv1 = sprintfc('%0.4e', errorsv1);
errors_scientificu2 = sprintfc('%0.4e', errorsu2);
errors_scientificv2 = sprintfc('%0.4e', errorsv2);
errors_scientificu3 = sprintfc('%0.4e', errorsu3);
errors_scientificv3 = sprintfc('%0.4e', errorsv3);

% % % % Display results in Tabular Format
table(M_values', errors_scientificu1', [ordersu1 NaN]', 'VariableNames', {'M', 'Erroru1', 'Order'})
table(M_values', errors_scientificv1', [ordersv1 NaN]', 'VariableNames', {'M', 'Errorv1', 'Order'})
table(M_values', errors_scientificu2', [ordersu2 NaN]', 'VariableNames', {'M', 'Erroru2', 'Order'})
table(M_values', errors_scientificv2', [ordersv2 NaN]', 'VariableNames', {'M', 'Errorv2', 'Order'})
table(M_values', errors_scientificu3', [ordersu3 NaN]', 'VariableNames', {'M', 'Erroru3', 'Order'})
table(M_values', errors_scientificv3', [ordersv3 NaN]', 'VariableNames', {'M', 'Errorv3', 'Order'})


