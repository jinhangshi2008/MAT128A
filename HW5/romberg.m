function intapprox =Jin_Hangshi_romberg(a,b,n,f)
%
% This function calculates the integral of f over [a,b] using Romberg
% Integration with degree n.
%
% Inputs:
%
% a, a real number
% b, a real number
% n, a positive integer
% f, a user-specified external function
%
% Outputs:
%
% intapprox, a real number
%
m=zeros(n,1);
for i=1:n
    m(i)=2.^(i-1);
end
h=(b-a)/m(1);
R=zeros(n,n);
R(1,1)=h*(f(a)+f(b))/2;
for i=2:n
    h=(b-a)/m(i);
    for j=2:m(i)
        R(i,1)=R(i,1)+h*(f(a+(j-1)*h));
    end
    R(i,1)=R(i,1)+h*(f(a)+f(b))/2;
end
for i=2:n
    for j=i:n
        R(j,i)=(4.^(i-1)*R(j,i-1)-R(j-1,i-1))/(4.^(i-1)-1);
    end
end
intapprox = vpa(R(n,n));
end