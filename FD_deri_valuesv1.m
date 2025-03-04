function[D1v, D2v, D3v] = FD_deri_valuesv1(x,v1)
N=length(x);
D1v=zeros(N,1); D2v=zeros(N,1); D3v=zeros(N,1);
LBC_support = 6; int_support = 4; RBC_support = 9;

for i=1:5

    % % % For First derivative Left Boundary approximation
    k=1; c1 = fdcoeffF(k,x(i),x(1:LBC_support));D1v(i)=c1*v1(1:LBC_support)';

    % % % For Second derivative Left Boundary approximation
    k=2; c2 = fdcoeffF(k,x(i),x(1:LBC_support));D2v(i)=c2*v1(1:LBC_support)';

    % % % For Third derivative Left Boundary approximation
    k=3; c3 = fdcoeffF(k,x(i),x(1:LBC_support));D3v(i)=c3*v1(1:LBC_support)';

end

for i= 6:N-6

    % % % For First derivative Interior approximation
    k=1; c1 = fdcoeffF(k,x(i),x(i-int_support:i+int_support));
         D1v(i)=c1*v1(i-int_support:i+int_support)';

    % % % For Second derivative Interior approximation
    k=2; c2 = fdcoeffF(k,x(i),x(i-int_support:i+int_support));
         D2v(i)=c2*v1(i-int_support:i+int_support)';

    % % % For Third derivative Interior approximation
    k=3; c3 = fdcoeffF(k,x(i),x(i-int_support:i+int_support));
         D3v(i)=c3*v1(i-int_support:i+int_support)';

end

for i=N-5:N

    % % % For First derivative Right Boundary approximation
    k=1; c1 = fdcoeffF(k,x(i),x(N-RBC_support:N));
         D1v(i)=c1*v1(N-RBC_support:N)';

    % % % For Seond derivative Right Boundary approximation
    k=2; c2 = fdcoeffF(k,x(i),x(N-RBC_support:N));
         D2v(i)=c2*v1(N-RBC_support:N)';

    % % % For Third derivative Right Boundary approximation
    k=3; c3 = fdcoeffF(k,x(i),x(N-RBC_support:N));
         D3v(i)=c3*v1(N-RBC_support:N)';

end

end