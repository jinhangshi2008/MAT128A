function px = newton(x,y)
n=size(x,1);                        %Get the number of rows in x
sum1=y(1);                          %Initiate sum with F0,0
b=sym('x');                         %Claim variable x
F=zeros(n,n);                       %Initiate matrix to store F
F(:,1)=y;                           %Let the first column be y
for i=2:n
    for j=2:i
        F(i,j)=(F(i,j-1)-F(i-1,j-1))/(x(i)-x(i-j+1));%Apply Algorithm 3.2
    end
    a=F(i,j);                       %Store the coef
    for k=2:i
    a=a*(b-x(k-1));                 %Multiply coef with each term
    end
    sum1=sum1+a;                    %Add up to P(x)
end
px= vpa(sum1,7);                    %Round coef to 7 digits
end
