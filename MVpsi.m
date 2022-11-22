function fval = MVpsi(l,xi,d,n)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
WM1 = zeros([1,d]);
WM2 = zeros([d,d]);
for i=1:n
    WM1 = WM1 + xi(i,:)*l(i)/n;
    WM2 = WM2 + xi(i,:).'*xi(i,:)*l(i)/n;
end
H = 2*WM2;
f = WM1;
Aeq = ones([1,d]);
beq = ones([1,1]);
A = [];
b = [];
lb = zeros([d,1]);
ub = [];
% x0 = ones([d,1])/d;
options = optimoptions(@quadprog, 'Display', 'none');
[x,fval] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
end