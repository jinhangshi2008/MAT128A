function px= LSQP(x,y,d)
A=zeros(d+1,d+1);                %Initiate matrix for coef of a
q=sym('x');                         %Claim variable x
for k=1:d+1
    b(k)=sum(y.*x.^(k-1));       %Calc y*x^(k-1)
    for i=k:d+k
        A(k,i-k+1)=sum(x.^(i-1));%Calc coef of a
    end
end
c=A\b';                  %Apply normal equation
s=0;
for j=1:d+1
    s=s+c(j)*q.^(j-1);
end
px = vpa(s,7);
end