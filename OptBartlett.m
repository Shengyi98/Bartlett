ns = [30,50];
qs = [1.642,2.7055,3.8415];
l_n = length(ns);
l_q = length(qs);
N=100000;
cover = zeros(l_n,l_q,N,3);
CIlength = zeros(l_n,l_q,N,3);

for idx_n = 1:l_n
    for idx_q = 1:l_q
        n = ns(idx_n);
        q0 = qs(idx_q); 
        truev = 1;
        
        % Estimate the theoretical A
        nt = 10000;
        Z = chi2rnd(2,1,nt);
        ep = randn(1,nt);
        Y = Z + ep;
        hatx = (Z*Z.')^(-1)*(Z*Y.');
        
        L = (Y-Z*hatx).^2;
        IF1 = L-mean(L);
        
        IF2 = zeros(nt);
        Elxx = 2*mean (Z.^2);
        
        for i=1:nt
            for j=1:nt
                IF2(i,j) = -8*Z(i)*Z(j)*(Y(i)-Z(i)*hatx)*(Y(j)-Z(j)*hatx)/Elxx;
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
        
        
        
        parfor rep = 1:N
        q = q0;
        
        Z = chi2rnd(4,1,n);
        ep = randn(1,n);
        Y = Z + ep;
        hatx = (Z*Z.')^(-1)*(Z*Y.');
        
        IF1 = zeros(n,1);
        L = (Y-Z*hatx).^2;
        IF1 = L-mean(L);
        
        IF2 = zeros(n);
        Elxx = 2*mean (Z.^2);
        
        for i=1:n
            for j=1:n
                IF2(i,j) = -4*Z(i)*Z(j)*(Y(i)-Z(i)*hatx)*(Y(j)-Z(j)*hatx)/Elxx;
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
        fun = @(l) (Y-Z*(Z.*l'*Z.')^(-1)*(Z.*l'*Y.')).*l'*(Y-Z*(Z.*l'*Z.')^(-1)*(Z.*l'*Y.')).'/n ;
        fun2 = @(l) -(Y-Z*(Z.*l'*Z.')^(-1)*(Z.*l'*Y.')).*l'*(Y-Z*(Z.*l'*Z.')^(-1)*(Z.*l'*Y.')).'/n ;
        l0 = ones(n,1);
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