ns = [15,30];
qs = [1.642,2.7055,3.8415];% (nominal level) quantile of chisquare 3.8415 (->95) 2.7055 (->90) 1.6424 (->80)
l_n = length(ns);
l_q = length(qs);
N=100000; % number of replications to estimate the coverage probability
shape = 2; % problem specific parameters
cover = zeros(l_n,l_q,N,3);%saves whether the true value is covered or not
CIlength = zeros(l_n,l_q,N,3);


h = @(x,y) min((x - y)^2 + x + y,12); % the objective function for V statistics

% estimate the true value of \psi(P) using sample average
truev = 0;
for i = 1:10000000
    r = randg(shape,2,1);
    truev = truev + (h(r(1),r(2))-truev)/i;
end

for idx_n = 1:l_n
    for idx_q = 1:l_q
        n = ns(idx_n);%number of samples per run
        q0 = qs(idx_q); 
        
        % estimate the coverage probability of different methods
        
        nt=5000;
        r = randg(shape,nt,1);
            
            H = zeros(nt);
            
            for i = 1:nt
                for j = 1:nt
                    H(i,j) = h(r(i),r(j));
                end
            end
            
            IF1 = zeros(nt,1);
            m = mean(H,'all');
            for i = 1:nt
                IF1(i) = 2*(mean(H(i,:))-m);
            end
            
            IF2 = zeros(nt);
            for i=1:nt
                for j=1:nt
                    %IF2(i,j) = H(i,j)-IF1(i)-IF1(j)-m;
                    IF2(i,j) = 2*(H(i,j)-IF1(i)/2-IF1(j)/2-m);
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
            % three parts of A(x)/-x where x^2=q
            %A1t = 1/9*q0^2*gamma^2/kappa^3;
            %A2t = 1/36*q0/kappa^3*(8*gamma^2-6*kappa*mu4+36*kappa^3+12*gamma*mu2c-12*gamma*kappa*mu2d);
            A3t = 1/36/kappa^3*(-12*gamma^2+18*kappa*mu4+36*kappa*(mu2a+2*mu2b)-36*gamma*mu2c-9*mu2c^2-18*kappa*mu2c*mu2d+9*kappa^2*(-2*mu22+mu2d^2-4*mu12d));
            At = A3t;
        
        
        parfor rep = 1:N
            
            
            
            r = randg(shape,n,1);
            
            H = zeros(n);
            
            for i = 1:n
                for j = 1:n
                    H(i,j) = h(r(i),r(j));
                end
            end
            
            IF1 = zeros(n,1);
            m = mean(H,'all');
            for i = 1:n
                IF1(i) = 2*(mean(H(i,:))-m);
            end
            
            IF2 = zeros(n);
            for i=1:n
                for j=1:n
                    IF2(i,j) = 2*H(i,j)-IF1(i)-IF1(j)-2*m;
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
            % three parts of A(x)/-x where x^2=q
            %A1 = 1/9*q0^2*gamma^2/kappa^3;
            %A2 = 1/36*q0/kappa^3*(8*gamma^2-6*kappa*mu4+36*kappa^3+12*gamma*mu2c-12*gamma*kappa*mu2d);
            A3 = 1/36/kappa^3*(-12*gamma^2+18*kappa*mu4+36*kappa*(mu2a+2*mu2b)-36*gamma*mu2c-9*mu2c^2-18*kappa*mu2c*mu2d+9*kappa^2*(-2*mu22+mu2d^2-4*mu12d));
            A = A3;
            
            
            
            fun = @(l) l.'*H*l/n/n ;
            fun2 = @(l) -l.'*H*l/n/n ;
            l0 = ones(n,1);
            % standard EL 
            q=q0;
            [~,fval_min_std] = fmincon(fun,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n));
            [~,fval_max_std] = fmincon(fun2,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n));
            % estimated Bartlett
            q = q0*(1+A/n);
            [~,fval_min_BL] = fmincon(fun,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n));
            [~,fval_max_BL] = fmincon(fun2,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n));
            % theoretical Bartlett
            q = q0*(1+At/n);
            [~,fval_min_TBL] = fmincon(fun,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n));
            [~,fval_max_TBL] = fmincon(fun2,l0,[],[],ones(1,n),n,zeros(n,1),[],@(l)KLdiv(l,q,n));
            
            % Record whether the true value is covered
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