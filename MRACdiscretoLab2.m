clear
clc
        num = 2.2439
        den = [0.1434 0.9837 1]
        theta = 0.19318
        ep = 1.2988;
        wn = 2.64074; 
        ts = (2*ep)/wn
        T = (0.2*(ts+theta)+0.6*(ts+theta))/2
        s = tf(num,den,'iodelay',theta)
        sd = c2d(s,T)
        [nd,dd] = tfdata(sd,'v')
        retd = totaldelay(sd)
        [z,p,b1] = tf2zpk(nd,dd)
       
        bneg = find(abs(z)>1)
        if isempty(bneg)
            bn = 1
        else
            bn = [1 -z(bneg)]
        end
        bpos = find(abs(z)<1)
        if isempty(bpos)
            bp = 1
        else
            bp = [1 -z(bpos)]
        end
        m = max(length(nd)-1,length(dd)-1);
        d = retd;
        if retd~=0 && mod(theta,T)~=0
            d = d-1;
            m = m+1;
        end
       

        zabs = exp(-ep*wn*T) % Valor absoluto de z
        theta2 = wn*T*sqrt(1-ep^2)
        [zx,zy] = pol2cart(theta2, zabs)
        ddm = conv([1 -(zx+zy*i)],[1 -(zx-zy*i)])
         
        syms z
        ams = vpa(subs(poly2sym(ddm(end:-1:1),z),z,z^-1),4)
        bnsim = vpa(subs(poly2sym(bn(end:-1:1),z),z,z^-1),4)
        k = sym2poly(vpa(1/limit(z^-1*bnsim/ams,z,1),4))
       
       
        % kd = expand(ams*z^2)
        % dd2 = coeffs(kd, z)
        syms z q0 q1 q2 p1
        p = '1';
        ordenpp = m+d-1-length(bpos); % Orden P'
        for i=1:1:ordenpp
            p = strcat(p,'+p',num2str(i),'*z^-',num2str(i));
        end
        p = str2sym(p)
        q = 'q0';
        ordenq = m-1;
        for j=1:1:ordenq
            q = strcat(q,'+q',num2str(j),'*z^-',num2str(j));
        end
        q = str2sym(q)
        A = vpa(subs(poly2sym(dd(end:-1:1),z),z^-1),4)
        B = vpa(subs(poly2sym(bn(end:-1:1),z),z^-1),4)
        ordenA = length(dd)-1      
        ordenBneg = length(bn)-1    
        n = max(ordenpp+ordenA,ordenq+ordenBneg+d+1)
        eq = vpa(expand(z^n*((p)*A+b1*(q)*B*z^-(d+1))),4)
        ecu = vpa(collect(eq),4)
        ecuv = vpa(coeffs(ecu,z),4)
        vec = zeros(1,n+1-length(ddm))
        sol = solve(ecuv(end:-1:1)==[ddm vec])
       
        solc = struct2cell(sol);
        psol(1) = 1;
        for i=1:ordenpp
            solc{i}
            psol(i+1)=sym2poly((vpa(solc{i},4)));
        end
        polinomioP = conv(bp, psol)
       
        for i=ordenpp+1:ordenpp+ordenq+1
            qsol(i-1)=sym2poly((vpa(solc{i},4)));
        end
       
        nd2 = k % Se cancelan las b menos en caso de hallar el modelo
        f = [(nd2*1)/(b1) 0] % El 1 es la referencia
       
        Ld = tf(nd,dd,T,'iodelay',d)
        Lr = tf(qsol,polinomioP,T)
        Slc = feedback(Ld,Lr)
       
        Lp = tf(f,polinomioP,T)
        Sfin = series(Slc,Lp)
        figure(5)
        step(Sfin)
