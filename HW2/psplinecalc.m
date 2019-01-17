function psplinecalc(x,y)
q=sym('x');                         %Claim variable x
n=size(x,1);                        %Get the number of rows in x
h=zeros(n-1,1);a=y;                 %Initiate vectors needed and store y into a
al=zeros(n-2,1);                    
l=zeros(n,1);u=zeros(n-1,1);
z=zeros(n,1);b=zeros(n-1,1);
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
    b(i)=(a(i+1)-a(i))/h(i)- h(i)*(c(i+1)+2*c(i))/3;
    d(i)=(c(i+1)-c(i))/(3*h(i));
end
a(n)=[];c(n)=[];                    %Remove the last entry of a and c
for i=1:n-1
    s{i}=vpa((a(i)+b(i)*(q-x(i))+c(i)*(q-x(i)).^2+d(i)*(q-x(i)).^3),7);
                                    %Store each natural cubic spline into s
end 
celldisp(s)                         %Display each natural cubic spline
end