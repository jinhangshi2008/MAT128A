function e = LSQPE(x,y,px)
n=size(x,1);                        %Get the number of rows in x
q=sym('x');                         %Claim variable x
s=0;
for i=1:n
    s=s+(y(i)-subs(px,q,x(i))).^2;  %Calc error
end
e=vpa(s,7);
end