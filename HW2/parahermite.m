function parahermite(x,y,a,b)
q=sym('t');                       %Claim variable x
x1=zeros(2,1);y1=zeros(2,1);      %Initiate vectors x1 and y1
x1(1)=a(1)-x(1);x1(2)=x(2)-a(2);  %Calculate and store alpha to x1 and beta to y1
y1(1)=b(1)-y(1);y1(2)=y(2)-b(2);
                                  %Compute x(t) and y(t) explicitly
xt=(2*(x(1)-x(2))+(x1(1)+x1(2)))*q.^3+(3*(x(2)-x(1))-(x1(2)+2*x1(1)))*q.^2+x1(1)*q+x(1)
yt=(2*(y(1)-y(2))+(y1(1)+y1(2)))*q.^3+(3*(y(2)-y(1))-(y1(2)+2*y1(1)))*q.^2+y1(1)*q+y(1)
end
