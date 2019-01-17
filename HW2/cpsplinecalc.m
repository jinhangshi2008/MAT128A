function cpsplinecalc(x,y,y1)
q=sym('x');                         %Claim variable x
n=size(x,1);                        %Get the number of rows in x
h=zeros(n-1,1);a=y;                 %Initiate vectors needed and store y into a
al=zeros(n+1,1);
l=zeros(n,1);u=zeros(n-1,1);
z=zeros(n,1);b=zeros(n-1,1);
c=zeros(n,1);d=zeros(n-1,1);
l(1)=1;FPO=y1(1);FPN=y1(2);         %Set the first entry of l to 1,
                                    %and set FPO to f'(x0), FPN to f'(xn)
for i=1:n-1
    h(i)=x(i+1)-x(i);               %Apply Algorithm 3.5
end
al(1)=3*(y(2)-y(1))/h(1)-3*FPO;     %Initiate the first and last entry
al(n)=3*FPN-3*(a(n)-a(n-1))/h(n-1); %of al
for i=2:n-1
    al(i)=3*(a(i+1)-a(i))/h(i)-3*(a(i)-a(i-1))/h(i-1);
end
l(1)=2*h(1);                        %Initiate the first entry of l,u, and z
u(1)=0.5;
z(1)=al(1)/l(1);
for i=2:n-1
    l(i)=2*(x(i+1)-x(i-1))-h(i-1)*u(i-1);
    u(i)=h(i)/l(i);
    z(i)=(al(i)-h(i-1)*z(i-1))/l(i);
end
l(n)=h(n-1)*(2-u(n-1));             %Set the last entry of l,z
z(n)=(al(n)-h(n-1)*z(n-1))/l(n);
c(n)=z(n);                          %store the last entry of z to the last of c
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