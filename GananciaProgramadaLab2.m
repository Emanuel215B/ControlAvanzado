
for var=1:1:5
    T=0.001;
    switch var
        case 1
            numerador=[0 1.6556];
            denominador=[0.42304 1];
            theta=0.38023;
        case 2
            numerador=[0 1.4644];
            denominador=[0.23512 1];
            theta=0.40501;
        case 3
            numerador=[0 2.2449];
            denominador=[0.84756 1];
            theta=0.35428;
        case 4
            numerador=[0 1.8225];
            denominador=[0.70652 1];
            theta=0.14038;
        case 5
            numerador=[0 0.84067];
            denominador=[1.3019 1];
            theta=0.32884;
    end
    sis=tf(numerador,denominador,'inputdelay',theta);
    sd=c2d(sis,T)
    [nd,dd]=tfdata(sd,'v');
    ret=totaldelay(sd);
    theta2=theta+T/2;
    ak=1.495;
    bk=-0.945;
    at=1.101;
    bt=0.771;
    ad=0.560;
    bd=1.006;
    k=numerador(1, 2);
    tao=denominador(1, 1);

    Kc=(ak/k)*(theta2/tao)^bk;
    ti=(tao/at)*(theta2/tao)^bt;
    td=(ad*tao)*(theta2/tao)^bd;
    q0(var)=Kc*(1+(T/(2*ti))+(td/T));
    q1(var)=-Kc*(1-(T/(2*ti))+((2*td)/T));
    q2(var)=(Kc*td)/T;
    nt(var,:)=nd;
    dt(var,:)=dd;
    retardos(var,1)=ret;

end



u=[38 42 46 50 54];
q0coe=polyfit(u,q0,4);
q1coe=polyfit(u,q1,4);
q2coe=polyfit(u,q2,4);

q0est=polyval(q0coe,u);
q1est=polyval(q1coe,u);
q2est=polyval(q2coe,u);

figure(1)
plot(u,q0est,u,q0,'+')

figure(2)
plot(u,q1est,u,q1,'+')
