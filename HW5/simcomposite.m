function px=simcomposite(a,b,n,f)
h=(b-a)/n;
sum=f(a)+f(b);
for i=1:n-2
    sum=sum+4*f(a+(i-1)*h)+2*f(a+i*h);
end
sum=sum+4*f(a+(n-1)*h);
px=sum*h/3;
end