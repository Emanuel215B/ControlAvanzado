clear;close;clc
[name,path]=uigetfile({'*.txt'},'Seleccione el archivo');
d=load([path name]);
d=d(1:end-1,:)
t=d(:,1);
y=d(:,2);
u=d(:,4);
figure(1)
plot(t,u,'r--',t,y,'b');
esc=46;
pos=find(u==esc);
ur=u(pos);
tr=t(pos);
yr=y(pos);
figure(2)
plot(tr,ur,'r--',tr,yr,'b--');
tt=tr-tr(1);
yt=yr-yr(1);
ut=ur-u(pos(1)-1);
figure(3)
plot(tt,ut,'r--',tt,yt,'b--');
%% MINIMOS CUADRADOS NO RECURSIVOS
U=ut';
Y=yt';
YT=-Y';
UT=U';
lu=length(U);
ly=length(Y);
n=2;
d=0;
l=n*2;
ceros=zeros(1,l);
F=ceros;

for i=1:1:ly
    F(i,1:end)=ceros;
end

for K=1:1:n
     F(K+1:end,K)=YT(1:ly-K);
end 

for K2=1:1:n
     F(K2+1+d:end,K2+n)=UT(1:ly-(K2+d));
end 

S=(inv(F'*F))*(F'*-YT);
num=zeros(1,n+1);
den=zeros(1,n+1);
num(1,1)=0;
den(1,1)=1;

for xnum=1:1:n
    num(1,xnum+1)=S(xnum+n,1);
    den(1,xnum+1)=S(xnum,1);
end

Fun=tf(num,den,1)
S2=step(UT(1)*Fun,1200);
figure(4)
plot(S2)
hold on
plot(Y)
hold on

%% MINIMOS CUADRADOS RECURSIVOS M1

numM1=zeros(1,n);
denM1=zeros(1,n);
numt=string(numM1);
dent=string(denM1);

for xnum=1:1:n
    x=S(xnum+n,1)+" Z^"+(-xnum);
    numt(1,xnum)=x;
end

for xden=1:1:n
    x2=S(xden,1)+" Z^"+(-xden);
    dent(1,1)="1";
    dent(1,xden+1)=x2;
end
%--------------------------------------------------------------------------
YN=8.7891;
% landa=input('Ingrese valor de lambda: ');
% a=input('Ingrese nuevo dato de a: ');

f2=zeros(2*n,1);

for h=1:1:n
    f2(h,1)=-Y(1,ly-h+1);
    f2(h+n,1)=U(1,lu-h+1);
end

PN=inv(F'*F);
LN=(1/1)*PN*f2*inv((1/1)+((1/1)*(f2'*PN*f2)));
SN=S+LN*(YN-f2'*S);

NS=zeros(1,n);
DS=zeros(1,n+1);
DS(1,1)=1;

for m=1:1:n
    NS(1,m)=SN(m+n,1);
    DS(1,m+1)=SN(m,1);
end 

Sol=tf(NS,DS,1)
SALIDA=step(UT(1)*Sol,1200);
figure(5)
plot(SALIDA,'B-')
hold on
plot(Y,'G')
hold on

%% MINIMOS CUADRADOS RECURSIVOS M2

theta=zeros(2*n,1);
I = eye(2*n);
alfa=1000;
P=alfa*I;
k1=n-1;
landa=1;
f=zeros(2*n,1);

for j=1:1:ly-2;
    for i=1:1:n
        f(i,1)=-Y(1,k1-i+2+(j-1));
        f(i+n,1)=U(1,k1-i+2-d+(j-1));
    end
    ft=f';
    e=(Y(k1+j+1)-ft*theta);
    L=(P*f)/(1+ft*P*f);
    theta2=theta+L*e;
    P2=1*(I-L*ft)*P;
    theta=theta2;
    P=P2;
end

numM2=zeros(1,n+1);
denM2=zeros(1,n+1);
numM2(1,1)=0;
denM2(1,1)=1;

for xnumM2=1:1:n
    numM2(1,xnumM2+1)=theta(xnumM2+n,1);
    denM2(1,xnumM2+1)=theta(xnumM2,1);
end

FunM2=tf(numM2,denM2,1)
S2M2=step(UT(1)*FunM2,1200);
figure(6)
plot(S2M2)
hold on
plot(Y,'g')
hold on