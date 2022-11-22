function [c,ceq] = KLdiv(l,q,n)


c1=0;
for i=1:n
    c1 = c1-2*log(l(i))+2*l(i)-2;
end
c = c1 - q;
ceq=[];
end