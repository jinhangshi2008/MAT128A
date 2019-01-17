function intvalue =Jin_Hangshi_splineint(x,y,a,b)
%
% This program calculates the integral of natural cubic splines of points (x,y)
% on [a,b].
%
% Inputs:
%
% x, a vector of length n+1
% y, a vector of length n+1
% a, a real number
% b, a real number
% 
% Outputs:
%
% intvalue, the value of the integral of a cubic spline
%
q=sym('x');                         %Claim variable x
n=size(x,1);                        %Get the number of rows in x
h=zeros(n-1,1);o=y;                 %Initiate vectors needed and store y into a
al=zeros(n-2,1);                    
l=zeros(n,1);u=zeros(n-1,1);
z=zeros(n,1);p=zeros(n-1,1);
c=zeros(n,1);d=zeros(n-1,1);
l(1)=1;                             %Set the first entry of l to 1
for i=1:n-1
    h(i)=x(i+1)-x(i);               %Apply Algorithm 3.4
end
for i=1:n-2
    al(i)=3*(y(i+2)-y(i+1))/h(i+1)-3*(y(i+1)-y(i))/h(i);
end
for i=2:n-1
    l(i)=2*(x(i+1)-x(i-1))-h(i-1)*u(i-1);
    u(i)=h(i)/l(i);
    z(i)=(al(i-1)-h(i-1)*z(i-1))/l(i);
end
l(n)=1;z(n)=0;                      %Set the last entry of l and 
                                    %z to 1 and 0 respectively
for i=n-1:-1:1
    c(i)=z(i)-u(i)*c(i+1);
    p(i)=(o(i+1)-o(i))/h(i)- h(i)*(c(i+1)+2*c(i))/3;
    d(i)=(c(i+1)-c(i))/(3*h(i));
end
o(n)=[];c(n)=[];                    %Remove the last entry of o and c
for i=1:n-1
    s{i}=(o(i)+p(i)*(q-x(i))+c(i)*(q-x(i)).^2+d(i)*(q-x(i)).^3);
                                    %Store each natural cubic spline into s
end                          
inte = int(s{1},[a,x(2)]);      %Initiate intvalue with the first int
                                    %from a to x(2)
for i=2:n-2
    inte = inte + int(s{i},[x(i),x(i+1)]);
                                    %Add up intvalue with more int from
                                    %x(i) to x(i+1)
end
intvalue = vpa(inte + int(s{n-1},[x(n-1),b]),7);
                                    %Add last int from x(n-1) to b
end