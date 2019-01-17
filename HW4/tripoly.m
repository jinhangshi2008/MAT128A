function sx=tripoly(fx,l,u,n)
a=zeros(n+1,1);b=zeros(n,1);
q=sym('x');
for i=1:n+1
    a(i)=int(fx*cos((i-1)*q),[l,u])/int(cos((i-1)*q).^2,[l,u]);
end
for i=2:n
    b(i)=int(fx*sin((i-1)*q),[l,u])/int(sin((i-1)*q).^2,[l,u]);
end
sum=0;
for i=2:n
    sum= sum+a(i)*cos((i-1)*q)+b(i)*sin((i-1)*q);
end
sx=vpa(a(1)/2+a(n+1)*cos(n*q)+ sum,7);
end