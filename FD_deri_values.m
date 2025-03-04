
function[u1,u2,u3] = FD_deri_values(x,u)
N=length(x);
u1=zeros(N,1); u2=zeros(N,1); u3=zeros(N,1);
LBC_support = 6; int_support = 4; RBC_support = 9;

for i=1:5

    % % % For First derivative Left Boundary approximation
    k=1; c1 = fdcoeffF(k,x(i),x(1:LBC_support));u1(i)=c1*u(1:LBC_support)';

    % % % For Second derivative Left Boundary approximation
    k=2; c2 = fdcoeffF(k,x(i),x(1:LBC_support));u2(i)=c2*u(1:LBC_support)';

    % % % For Third derivative Left Boundary approximation
    k=3; c3 = fdcoeffF(k,x(i),x(1:LBC_support));u3(i)=c3*u(1:LBC_support)';

end

for i= 6:N-6

    % % % For First derivative Interior approximation
    k=1; c1 = fdcoeffF(k,x(i),x(i-int_support:i+int_support));
         u1(i)=c1*u(i-int_support:i+int_support)';

    % % % For Second derivative Interior approximation
    k=2; c2 = fdcoeffF(k,x(i),x(i-int_support:i+int_support));
         u2(i)=c2*u(i-int_support:i+int_support)';

    % % % For Third derivative Interior approximation
    k=3; c3 = fdcoeffF(k,x(i),x(i-int_support:i+int_support));
         u3(i)=c3*u(i-int_support:i+int_support)';

end

for i=N-5:N

    % % % For First derivative Right Boundary approximation
    k=1; c1 = fdcoeffF(k,x(i),x(N-RBC_support:N));
         u1(i)=c1*u(N-RBC_support:N)';

    % % % For Seond derivative Right Boundary approximation
    k=2; c2 = fdcoeffF(k,x(i),x(N-RBC_support:N));
         u2(i)=c2*u(N-RBC_support:N)';

    % % % For Third derivative Right Boundary approximation
    k=3; c3 = fdcoeffF(k,x(i),x(N-RBC_support:N));
         u3(i)=c3*u(N-RBC_support:N)';

end

end


