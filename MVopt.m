ns = [15,30,50];
qs = [1.642,2.7055,3.8415];
l_n = length(ns);
l_q = length(qs);
N = 5000;
d = 3;
cover = zeros(l_n,l_q,N,3);
CIlength = zeros(l_n,l_q,N,3);

truev = 1/d;
for idx_n = 1:l_n
    for idx_q = 1:l_q
        n = ns(idx_n);
        q0 = qs(idx_q); 
        % Estimate the theoretical A
        nt = 1000; % sample size; needs to be large
        xi = normrnd(0,1,[nt,d]);
        %ep = randn(1,nt);
        %Y = Z + ep;
        %hatx = (Z*Z.')^(-1)*(Z*Y.');
        
        IF1 = zeros([nt,1]);
        for i=1:nt
            me = mean(xi(i,:));
            IF1(i)= me + me^2 - 1/d;
        end
        
        IF2 = zeros(nt);
        
        for i=1:nt
            for j=1:nt
                me1 = mean(xi(i,:));
                me2 = mean(xi(j,:));
                xi12 = dot(xi(i,:),xi(j,:));
                IF2(i,j) = -(1+2*me1)*(1+2*me2)*(xi12-d*me1*me2);
            end
        end
        
        kappa = mean(IF1.^2);
        gamma = mean(IF1.^3);
        mu4 = mean(IF1.^4);
        
        K = zeros(nt,1);
        for i =1:nt
            for j = 1:nt
                K(i) = K(i)+IF2(i,j)*IF1(j);
            end
            K(i) = K(i)/nt;
        end
        mu2a = mean(K.^2);
        mu2b = 0;
        mu2c = 0;
        
        for i = 1:nt
            for j = 1:nt
                mu2b = mu2b+IF1(i)*IF1(j)^2*IF2(i,j);
                mu2c = mu2c+IF1(i)*IF1(j)*IF2(i,j);
            end
        end
        
        mu2b = mu2b/nt/nt;
        mu2c = mu2c/nt/nt;
        mu22 = mean(IF2.^2,'all');
        mu12d = 0;
        mu2d = 0;
        for i = 1:nt
            mu12d = mu12d + IF1(i)*IF2(i,i);
            mu2d = mu2d + IF2(i,i);
        end
        
        mu12d = mu12d/nt;
        mu2d = mu2d/nt;
        
        At = 1/36/kappa^3*(-12*gamma^2+18*kappa*mu4+36*kappa*(mu2a+2*mu2b)-36*gamma*mu2c-9*mu2c^2-18*kappa*mu2c*mu2d+9*kappa^2*(-2*mu22+mu2d^2-4*mu12d));
        
        
        tic
        parfor rep = 1:N
        
        xi = normrnd(0,1,[n,d]);
        l0 = ones(n,1);
        xhat = MVpsiarg(l0,xi,d,n);
        IF1 = zeros(n,1);
        for i=1:n
            mehat = dot(xi(i,:),xhat);
            IF1(i) = mehat + mehat^2;
        end
        IF1 = IF1 - mean(IF1);
        
        IF2 = zeros(n);
        Elxx = zeros(d-1);
        for i=1:n
            veci = xi(i,1:(d-1))-xi(i,d)*ones(1,d-1);
            Elxx = Elxx + veci.'*veci;
        end
        Elxxinv = inv(Elxx)*n;
        
        for i=1:n
            for j=1:n
                vec1 = xi(i,1:(d-1))-xi(i,d)*ones(1,d-1);
                vec2 = xi(j,1:(d-1))-xi(j,d)*ones(1,d-1);
                lx1 = (1+2*(dot(vec1,xhat(1:(d-1)))+xi(i,d)))*vec1;
                lx2 = (1+2*(dot(vec2,xhat(1:(d-1)))+xi(j,d)))*vec2;
                IF2(i,j) = -2*lx1*Elxxinv*lx2.';
            end
        end
        
        kappa = mean(IF1.^2);
        gamma = mean(IF1.^3);
        mu4 = mean(IF1.^4);
        
        K = zeros(n,1);
        for i =1:n
            for j = 1:n
                K(i) = K(i)+IF2(i,j)*IF1(j);
            end
            K(i) = K(i)/n;
        end
        mu2a = mean(K.^2);
        mu2b = 0;
        mu2c = 0;
        
        for i = 1:n
            for j = 1:n
                mu2b = mu2b+IF1(i)*IF1(j)^2*IF2(i,j);
                mu2c = mu2c+IF1(i)*IF1(j)*IF2(i,j);
            end
        end
        
        mu2b = mu2b/n/n;
        mu2c = mu2c/n/n;
        mu22 = mean(IF2.^2,'all');
        mu12d = 0;
        mu2d = 0;
        for i = 1:n
            mu12d = mu12d + IF1(i)*IF2(i,i);
            mu2d = mu2d + IF2(i,i);
        end
        
        mu12d = mu12d/n;
        mu2d = mu2d/n;
        
        A = 1/36/kappa^3*(-12*gamma^2+18*kappa*mu4+36*kappa*(mu2a+2*mu2b)-36*gamma*mu2c-9*mu2c^2-18*kappa*mu2c*mu2d+9*kappa^2*(-2*mu22+mu2d^2-4*mu12d));
        A ;
        
        
        options = optimoptions('fmincon','MaxFunctionEvaluations',10000);
        fun = @(l) MVpsi(l,xi,d,n);
        fun2 = @(l) -MVpsi(l,xi,d,n) ;
        
        q=q0;
        [x_min,fval_min_std] = fmincon(fun,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n),options);
        [x_max,fval_max_std] = fmincon(fun2,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n),options);
        
        q = q0*(1+A/n);
        [x_min,fval_min_BL] = fmincon(fun,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n),options);
        [x_max,fval_max_BL] = fmincon(fun2,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n),options);
        
        q = q0*(1+At/n);
        [x_min,fval_min_TBL] = fmincon(fun,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n),options);
        [x_max,fval_max_TBL] = fmincon(fun2,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n),options);
        
        cover1 = 0;
        cover2 = 0;
        cover3 = 0;
        
        length1 = -fval_max_std - fval_min_std;
        length2 = -fval_max_BL - fval_min_BL;
        length3 = -fval_max_TBL - fval_min_TBL;
        
        if fval_min_std<truev && truev< -fval_max_std
            cover1 = 1;
        end
        
        if fval_min_BL<truev && truev< -fval_max_BL
            cover2 = 1;
        end
        
        if fval_min_TBL<truev && truev< -fval_max_TBL
            cover3 = 1;
        end
        
        cover(idx_n,idx_q,rep,:) = [cover1,cover2,cover3];
        CIlength(idx_n,idx_q,rep,:) = [length1,length2,length3];
        [rep,A,At]
        end
    end
end
save('cover','cover');
save('length','CIlength');