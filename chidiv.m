function [c,ceq] = chidiv(l,q,n)

c1=0;
for i=1:n
    c1 = c1+(l(i)-1)^2;
end
c = c1 - q;
ceq=[];
end