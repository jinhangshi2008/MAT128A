function px=composite(a,b,n,f)
h=(b-a)/n;       %Calculate step size
sum=f(a)+f(b);   %Initiate sum with f(a)+f(b)
for i=2:n
    sum=sum+2*f(a+(i-1)*h);
                 %Apply composite trapzoid rule
end
px=sum*h/2;      %Return the result
end