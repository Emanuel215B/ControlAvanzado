num=[0.48828 0.043592];
den=[1 -0.71672 -0.24061];

m=length(num)-1;
d=1;
u=m+d;
v=m;

Funcion=tf(num,den,1,'var','z^-1');

l=max([m+u,m+d+v]);
matriz=2*m+d+1;

M = zeros(matriz);

numt=num';
dent=den';

dent(1)=[];
la=length(dent); 

numt(1)=[];
lb=length(numt);

for q=1:1:m+d
    M(q,q)=1;
    M(q+1:la+q,q)=dent;
    M(matriz,q)=1;
end 

for q2=1:1:m+1
    M(d+q2:lb+d+q2-1,m+d+q2)=numt;
end 

alpha=poly([0.7+0.3i,0.7-0.3i])
alpha(1)=[];
lalph=length(alpha);
Vector=zeros(2*m+d+1,1);

Vector(2*m+d+1,1)=-1;
for q3=1:1:lalph
    Vector(q3,1)=alpha(q3)-dent(q3);
end 

PQ=inv(M)*Vector;

numDZT=PQ(m+d+1:end);
denDZT=PQ(1:m+d);
numDZ=numDZT';
denDZ=[1;denDZT]';

DZ=tf(numDZ,denDZ,1,'var','z^-1')

ST=series(DZ,Funcion);
SF=feedback(ST,tf(1,1));
step(SF*4);

