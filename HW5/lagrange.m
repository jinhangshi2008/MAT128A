function pc = lagrange(x,y,c)
% This function is to calculate the P(x) using
% lagrange interpolation at x=c.
% 
% Inputs:
%
% x, a vector of length n+1
% y, a vector of length n+1
% c, a real number
%
% Outputs:
%
% pc, the value of P(c)
n=size(x,1);                        %Get the number of rows in x
sum1=0;                             %Initiate sum with 0
b=sym('x');                         %Claim variable x
for j=1:n
    a=x;
    a(j)=[];                        %For each iteration, exclude the jth entry
    L=expand(prod((b-a)./(x(j)-a)));%Get each L(x) and expand the polynomial
    sum1=sum1+L*y(j);               %Add each L(x) to get P(x)
end
pc = vpa(subs(sum1,b,c));           %Plug in c to get P(c).
end