% % % % "RBF based backward differentiation methods for stiff differential equations"
% % % % Authors: A. Sreedhar, Manoj Kumar Yadav, Chirala Satyanarayana


% % % % To solve the Non linear system of stiff equations using BDF methods
% % % % Define the Sytem of nonlinear ODEs of the form du/dt= −1002 u(t) + 1000 v(t)^2 and 
% % % %                                                dv/dt= u(t) − v(t)(1 + v(t))
% % % % w.r.t. initial conditions u(0)=1; and v(0)=1;

clc
clear all

% % % % % % Parameters
a = 0; b = 5; M = 51;  % % % %  Number of time steps
t = linspace(a, b, M); % % % %  Time vector
u0 = 1; v0 = 1; % % % %  Initial conditions
dt = (b-a)/(M-1); % % % % Step size
% % % % % % Define the exact solutions
uexact = @(t) exp(-2*t);% % % % Exact Solution for u
vexact = @(t) exp(-t);  % % % % Exact Solution for v

% % % % Setting Initial conditions
u(1) = u0; v(1) = v0;

% % % % Defining the Given vector for RHS function
F = @(u, v) [-1002*u+1000*v^2; u-v*(1+v)];

% % % % Estimating the u(2) solution
u(2) = uexact(t(2)); v(2) = vexact(t(2));

% % % % Now apply BDF2 for subsequent steps
for i = 3:M
    u(i) = u(i-1); % % % % Initial guess for u(n+1)
    v(i) = v(i-1); % % % % Initial guess for v(n+1)

    % % % % Newton's method to solve the nonlinear system for u(n+1) and v(n+1)
    tol2=exp(-16); err2=10;
    while err2>tol2 
        FC1 = ((3/2)*u(i)-2*u(i-1)+(1/2)*u(i-2))/(dt)+1002*u(i)-1000*v(i)^2;
        FC2 = ((3/2)*v(i)-2*v(i-1)+(1/2)*v(i-2))/(dt)-u(i)+v(i)*(1+v(i));

        JC11 = (3/(2*dt))+1002; JC12 = -2000*v(i); 
        JC21 = -1; JC22 = (3/(2*dt))+(1+2*v(i));
        J1 = [JC11, JC12; JC21, JC22]; % % % % Jacobian matrix

        Fvec2 = [FC1; FC2];
        delta2 = -J1\ Fvec2;  % % % % Newton's method update
        u(i) = u(i)+delta2(1);
        v(i) = v(i)+delta2(2);
        err2=norm(Fvec2,inf);
    end  
end 

% % % % % Call the derivative values by Fornberg algorithm
[D1u, D2u, D3u] = FD_deri_valuesu1(t, u);
[D1v, D2v, D3v] = FD_deri_valuesv1(t, v);

% % % % % For choosing the optimised Epsilon Value Intialisation
epg1=zeros(1,M); epg2=zeros(1,M);
epm1=zeros(1,M); epm2=zeros(1,M);

% % % % % To find the u(2) and v(2) solution using Classical BDF1 Scheme
u1(1) = u0;   v1(1) = v0; 
u1(2) = u(2); v1(2) = v(2);

% % % % % Now apply RBF-BDF2 for subsequent steps
for i = 3:M
    u1(i) = u1(i-1); % Initial guess for u(n+1)
    v1(i) = v1(i-1); % Initial guess for v(n+1)

    % Newton's method to solve the nonlinear system for u(n+1) and v(n+1)
    tol=exp(-16); err=10;
    while err>tol      

%         % % % % % Calculate epsilon(i) for GA-BDF2
%         epg1(i) = -(D3u(i))/(6*D1u(i));  % % opt. epsilon with GA-Weights
%         epg2(i) = -(D3v(i))/(6*D1v(i));  % % opt. epsilon with GA-Weights
%         
%         % % % % % % Define GA-BDF2 u_weights
%         alphag22 = (3/(2*dt))-(11*epg1(i)*dt)/6+(37*epg1(i)^2*dt^3)/36+(445*epg1(i)^3*dt^5)/216-(2129*epg1(i)^4*dt^7)/1296;
%         alphag21 = (-2/dt)+(5*epg1(i)*dt)/3-(13*epg1(i)^2*dt^3)/18-(445*epg1(i)^3*dt^5)/108-(3631*epg1(i)^4*dt^7)/648;
%         alphag20 = (1/(2*dt))+(epg1(i)*dt)/6-(11*epg1(i)^2*dt^3)/36+(445*epg1(i)^3*dt^5)/216+(9391*epg1(i)^4*dt^7)/1296;
%         
%         % % % % Define GA-BDF2 v_weights
%         betag22 = (3/(2*dt))-(11*epg2(i)*dt)/6+(37*epg2(i)^2*dt^3)/36+(445*epg2(i)^3*dt^5)/216-(2129*epg2(i)^4*dt^7)/1296;
%         betag21 = (-2/dt)+(5*epg2(i)*dt)/3-(13*epg2(i)^2*dt^3)/18-(445*epg2(i)^3*dt^5)/108-(3631*epg2(i)^4*dt^7)/648;
%         betag20 = (1/(2*dt))+(epg2(i)*dt)/6-(11*epg2(i)^2*dt^3)/36+(445*epg2(i)^3*dt^5)/216+(9391*epg2(i)^4*dt^7)/1296;
%         
%         F1 = alphag22*u1(i)+alphag21*u1(i-1)+alphag20*u1(i-2)+1002*u1(i)-1000*v1(i)^2;
%         F2 = betag22*v1(i)+betag21*v1(i-1)+betag20*v1(i-2)-u1(i)+v1(i)*(1+v1(i));
%         
%         Jg11 = alphag22+1002; Jg12 = -2000*v1(i);
%         Jg21 = -1; Jg22 = betag22+(1+2*v1(i));
%         J = [Jg11, Jg12; Jg21, Jg22];  % Jacobian matrix

        % % % % Calculate epsilon(i) for MQ-BDF2
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

        F1 = alpham22*u1(i)+alpham21*u1(i-1)+alpham20*u1(i-2)+1002*u1(i)-1000*v1(i)^2;
        F2 = betam22*v1(i)+betam21*v1(i-1)+betam20*v1(i-2)-u1(i)+v1(i)*(1+v1(i));

        Jm11 = alpham22+1002; Jm12 = -2000*v1(i);
        Jm21 = -1; Jm22 = betam22+(1+2*v1(i));
        J = [Jm11, Jm12; Jm21, Jm22]; % % % % % Jacobian matrix

        Fvec = [F1; F2];
        delta = -J\ Fvec; % % % % % Newton's method update
        u1(i) = u1(i)+delta(1);
        v1(i) = v1(i)+delta(2);
        err=norm(Fvec,inf);
    end
end

% % % % % To find the errors with classical BDF Schemes for u and v
erroruc2 = abs(uexact(t(M))-u(M))
errorvc2 = abs(vexact(t(M))-v(M))

% % % % % To find the errors with RBF-BDF Schemes for u and v
erroru2 = abs(uexact(t(M))-u1(M))
errorv2 = abs(vexact(t(M))-v1(M))
