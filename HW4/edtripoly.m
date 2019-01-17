function E=edtripoly(fx,l,u,n,m)
a=zeros(n+1,1);b=zeros(n,1);       %Initiate a_k and b_k
q=sym('x');                        %Claim variable x
x=zeros(2*m,1);y=zeros(2*m,1);     %Initiate x_j and y_j
for i=1:2*m
    x(i)=l+((i-1)/m)*u;            %Calculate x_j
end
for i=1:2*m
    y(i)=subs(fx,q,x(i));          %Calculate y_j
end
for i=1:n+1
    for j=1:2*m
        a(i)=a(i)+y(j)*cos((i-1)*x(j));
    end
    a(i)=a(i)/m;                   %Calculate a_k
end
for i=2:n
    for j=1:2*m
    b(i)=b(i)+y(j)*sin((i-1)*x(j));
    end
    b(i)=b(i)/m;                   %Calculate b_k
end
sum=0;                             %Initiate sum
for i=2:n
    sum= sum+a(i)*cos((i-1)*q)+b(i)*sin((i-1)*q);
                                   %Calculate the sum for k=1,...,n-1
end
sx=a(1)/2+a(n+1)*cos(n*q)+ sum;
                                   %Calculate S_n(x)
E=0;                               %Initiate error
for i=1:2*m
    E=E+(y(i)-subs(sx,q,x(i))).^2; %Calculate error
end
E=vpa(E,7);
end