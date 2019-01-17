function hx = Hermite(x,y,y1)
b=sym('x');                         %Claim variable x
n=size(x,1);                        %Get the number of rows in x
X=[x(1) x(1)];                      %Initiate the first entries of vectors 
Y=[y(1) y(1)];                      %X and Y with the first entries of x and y
for i=2:n
    X=[X x(i) x(i)];                %Complete X and Y by copying each entry 
    Y=[Y y(i) y(i)];                %and inserting it right below the one copied
end
sum1=y(1)+y1(1)*(b-x(1));           %Initiate sum at F(0,0)
F=zeros(2*n,2*n);                   %Initiate matrix to store F
F(:,1)=Y';                          %Store the first column of F with Y
Y1=[0 y1(1)];                       %Initiate Y1 with first entry of y1
for i=2:n
    Y1=[Y1 (y(i)-y(i-1))/(x(i)-x(i-1)) y1(i)];%Construct Y1 with each of y1 on
                                              %odd position and the
                                              %computation on even.
end
F(:,2)=Y1';                         %Store the second column of F with Y1
for i=3:2*n
    for j=3:i
        F(i,j)=(F(i,j-1)-F(i-1,j-1))/(X(i)-X(i-j+1));%Apply Algorithm 3.3
    end
    a=F(i,j);                       %Store the coef
    for k=2:i
    a=a*(b-X(k-1));                 %Multiply coef with each term
    end
    sum1=sum1+a;                    %Add up to H(x)
end
hx= vpa(sum1,7);                    %Round coef to 7 digits
end