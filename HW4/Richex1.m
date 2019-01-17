function N=Richex1(fx,d,h,x0)
N=zeros(d,d);                   %Initiate N as a dXd matrix
q=sym('x');                     %Claim variable x
for i=1:d
    N(i,1)=(subs(fx,q,x0+(h/(2.^(i-1))))-subs(fx,q,x0))/(h/(2.^(i-1)));
                                %Calculate O(h)
end
for j=2:d
    for i=j:d
         N(i,j)=(1/(2.^(j-1)-1))*(2.^(j-1)*N(i,j-1)-N(i-1,j-1));
                                %Calculate O(h^j)
    end
end
end
