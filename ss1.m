function f=ss(x,p)


alpha=p(1);
beta=p(2);
gamma=p(3);
a=p(4);
lambda=p(5);
theta=p(6);
niu=p(7);
b=p(8);
d=p(9);
sigmac=p(10);
sigma=p(11);
omega=p(12);
tau=p(13);
phi=p(14);
betadisc=p(15);
delta=p(16);
nu=p(17);
poaverage=p(18);
A=p(19);
rho=p(20);
RBAR=p(21);
bss=p(22);
psi=p(23);
gexg=p(24);
oexg=p(25);
igexg=p(26);
posubexg=p(27);
integcost=p(28);

%variable definitions

OE     	=	 x(1);
E    	=	 x(2);
G     	=	 x(3);
RNW   	=	 x(4);
PE    	=	 x(5);
POSUB 	=	 x(6);
PG    	=	 x(7);
PRNW  	=	 x(8);
S     	=	 x(9);
OS    	=	 x(10);
PS     	=	 x(11);
N     	=	 x(12);
KP     	=	 x(13);
SF   	=	 x(14);
W     	=	 x(15);
R    	=	 x(16);
SH    	=	 x(17);
C     	=	 x(18);
RSTAR 	=	 x(19);
TR    	=	 x(20);
PO    	=	 x(21);
B     	=	 x(22);
O     	=	 x(23);
KG    	=	 x(24);
Y    	=	 x(25);




f(1)=E-alpha*OE-beta*G-gamma*RNW;
f(2)=PE-1/alpha*POSUB;
f(3)=PG-beta/alpha*POSUB;
f(4)=PRNW-gamma/alpha*POSUB;
f(5)=S-(a*E^lambda+(1-a)*OS^lambda)^(1/lambda);
f(6)=PS-POSUB*(OS^(1-lambda))/(1-a)*(a*E^lambda+(1-a)*OS^lambda)^((lambda-1)/lambda);
f(7)=PS-PE*(E^(1-lambda))/a*(a*E^lambda+(1-a)*OS^lambda)^((lambda-1)/lambda);
f(8)=W-theta*N^(theta-1)*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu);
f(9)=R-(1-theta)*(1-b)*KP^(niu-1)*N^theta*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu-1);
f(10)=PS-(1-theta)*b*SF^(niu-1)*N^theta*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu-1);
f(11)=N-1;%W-omega*N^(1/tau)/(d*SH^(sigmac-1)*(C^sigmac+d*SH^sigmac)^((1-sigmac-sigma)/sigmac));
f(12)=PS-d*(SH/C)^(sigmac-1);
f(13)=(B-(1+RSTAR)*B)-bss*(Y+PO*(O-OE-OS));
f(14)=1-betadisc*(1-delta+R);
f(15)=B+KP-(1-delta)*KP+C+PS*SH-W*N-R*KP-(1+RSTAR)*B-TR;
f(16)=PRNW*RNW+POSUB*(OS+OE)+PO*(O-OS-OE)+PG*G-TR-igexg;
f(17)=RNW-A/(1+integcost*RNW/E)*KG;
f(18)=S-SF-SH;
f(19)=PO-poaverage;
f(20)=RSTAR-RBAR;
f(21)=nu*KG-igexg;
f(22)=POSUB-posubexg;
f(23)=O-oexg;
f(24)=G-gexg;
f(25)=Y-N^theta*((1-b)*KP^niu+b*SF^niu)^((1-theta)/niu);


f=f';