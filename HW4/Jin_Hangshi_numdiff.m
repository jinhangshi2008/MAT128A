function derivs = Jin_Hangshi_numdiff(a,h,f)
%
% This function calculates f'(a) by using first the three-point midpoint and
% endpoint formulas with step size h, then the three-point endpoint formula with
% step size -h, and finally the five-point midpoint formula with step size h.
% 
% Inputs:
%
% a, a real number
% h, a positive real number
% f, a user-specified external function
%
% Outputs:
%
% derivs, a vector of derivative approximations
%
df=zeros(size(a,2),4);              
df(:,1)=(f(a+h)-f(a-h))/(2*h);
df(:,2)=(-3*f(a)+4*f(a+h)-f(a+2*h))/(2*h);
df(:,3)=(-3*f(a)+4*f(a-h)-f(a-2*h))/(-2*h);
df(:,4)=(f(a-2*h)-8*f(a-h)+8*f(a+h)-f(a+2*h))/(12*h);
derivs = df;
end