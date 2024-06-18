%% LABELING BLOCK
var c, l, n, nf, np, kf, kp, i, if, ip, y, yf, yp, v, A, % endogenous variables (quantities)
    lamf, lamp, gamf, gamp, % endogenous variables (prices)
    zf, zp, pk, pv, r, % exogenous variables
    nx, Lp; % additional calculations

varexo ezf, ezp, epk, epv, er;

%% PARAMETER BLOCK

parameters BETA, UPSILON, SIGMA, PSI, %prefrences
           THETA1, THETA2, ZETAP, %technology of primary goods sector
           ETAP, ETAF,
           S, SK, SV,
           ALPHA, U, DELTA, TB, ZETAF,
           rhozf, rhozp, rhopk, rhopv, rhor;

BETA=0.97;
UPSILON=(1+0.83)/0.83;
SIGMA=2.61;
PSI=5.35;

THETA1=0.37;
THETA2=0.18;
ETAP=2.2;
ZETAP=(ETAP-1)/ETAP;

ALPHA=0.45;

S=0.55;
SK=0.23;
SV=0.32;

U=(1-0.77)/0.77;

DELTA=0.10;
TB=-0.096;
ETAF=2;
ZETAF=(ETAF-1)/ETAF;

rhopk=0.44;
rhopv=0.42;
rhor=0.34;
rhozf=0.54;
rhozp=0.32;

%% MODEL BLOCK
model; 

% Model equations
(c-PSI*(1-l)^UPSILON)^(-SIGMA)=lamf; 
(PSI*UPSILON*(1-l)^(UPSILON-1))*((c-PSI*(1-l)^UPSILON)^(-SIGMA))=lamf*(ALPHA)*zf*(nf)^(ALPHA-1)*(S*kf^(-U)+(1-S)*v^(-U))^(-(1-ALPHA)/U);
-gamf*DELTA^(1-ZETAF)*(if/kf)^(ZETAF-1)*kf=lamf;
-gamp*DELTA^(1-ZETAP)*(ip/kp)^(ZETAP-1)*kp=lamf*pk;
lamf*(1-ALPHA)*(1-S)*zf*nf^(ALPHA)*(S*kf^(-U)+(1-S)*v^(-U))^((-(1-ALPHA)/U)-1)*v^(-U-1)=lamp*pv;
lamf*ALPHA*zf*nf^(ALPHA-1)*(S*kf^(-U)+(1-S)*v^(-U))^(-(1-ALPHA)/U)=lamp*THETA1*zp*np^(THETA1-1)*kp^(THETA2)*Lp^(1-THETA1-THETA2);
gamf(-1)=BETA*lamf*S*(1-ALPHA)*zf*nf^(ALPHA)*(S*kf^(-U)+(1-S)*v^(-U))^((-(1-ALPHA)/U)-1)*kf^(-U-1)+BETA*gamf*((1-DELTA)+DELTA^(1-ZETAF)*(if/kf)^(ZETAF-1)*(if/kf)+(1/ZETAF)*(DELTA^(1-ZETAF)*(if/kf)^(ZETAF)-(1-ZETAF)*DELTA));
gamp(-1)=BETA*lamp*THETA2*zp*np^(THETA1)*kp^(THETA2-1)*Lp^(1-THETA1-THETA2)+BETA*gamp*((1-DELTA)+DELTA^(1-ZETAP)*(ip/kp)^(ZETAP-1)*(ip/kp)+(1/ZETAP)*(DELTA^(1-ZETAP)*(ip/kp)^(ZETAP)-(1-ZETAP)*DELTA));
lamp(-1)=BETA*lamp*(1+r);

%Additional Equations
yf=zf*nf^(ALPHA)*(S*kf^(-U)+(1-S)*v^(-U))^(-(1-ALPHA)/U);
yp=zp*np^(THETA1)*kp^(THETA2)*Lp^(1-THETA1-THETA2);
kf=(1-DELTA)*kf(-1)+(1/ZETAF)*(DELTA^(1-ZETAF)*(if(-1)/kf(-1))^(ZETAF)-(1-ZETAF)*DELTA)*kf(-1);
kp=(1-DELTA)*kp(-1)+(1/ZETAP)*(DELTA^(1-ZETAP)*(ip(-1)/kp(-1))^(ZETAP)-(1-ZETAP)*DELTA)*kp(-1);
c+if=yf;
nx=yp-pk*ip-pv*v;
l+nf+np=1;
A=nx(-1)+A(-1)*(1+r(-1));
y=yf+yp;
n=nf+np;
i=if+ip;

% Exogenous processes
log(zf)=rhozf*log(zf(-1))+ezf;
log(zp)=rhozp*log(zp(-1))+ezp;
log(pk)=rhopk*log(pk(-1))+epk;
log(pv)=rhopv*log(pv(-1))+epv;
r=rhor*r(-1)+er;

end;

%%Steady State block
zfstar=1;
zpstar=1;
pkstar=1;
pvstar=1;
rstar=(1/BETA)-1;
Astar=TB;
nxstar=-rstar*Astar;
ypstar=(BETA*(1+DELTA)-1)/(BETA*THETA2);
kpstar=THETA2*ypstar;
npstar=THETA1*ypstar;
ipstar=DELTA*kpstar;
vstar=(1-DELTA*THETA2)*ypstar-nxstar;
yfstar=(1/SV)*vstar;
nfstar=ALPHA*yfstar;
kfstar=SK*yfstar;
ifstar=DELTA*kfstar;
lstar=1-nfstar-npstar;
Lpstar=((npstar)^(-THETA1)*(kpstar)^(-THETA2)*ypstar)^(1/(1-THETA1-THETA2));
cstar=(1-DELTA*SK)*yfstar;
lamfstar=(cstar-PSI*(1-lstar)^UPSILON)^(-SIGMA);
lampstar=lamfstar*(1-ALPHA)*(1-S)*nfstar^(ALPHA)*(S*kfstar^(-U)+(1-S)*vstar^(-U))^((-(1-ALPHA)/U)-1)*vstar^(-U-1);
gamfstar=-lamfstar*kfstar^(-1);
gampstar=-lampstar*kpstar^(-1);
ystar=yfstar+ypstar;
nstar=nfstar+npstar;
istar=ifstar+ipstar;

initval;
c=cstar;
l=lstar;
nf=nfstar;
np=npstar;
kf=kfstar;
kp=kpstar;
if=ifstar;
ip=ipstar;
yf=yfstar;
yp=ypstar;
v=vstar;
A=Astar;
lamf=lamfstar;
lamp=lampstar;
gamf=gamfstar;
gamp=gampstar;
zf=zfstar;
zp=zpstar;
pk=pkstar;
pv=pvstar;
r=rstar;
nx=nxstar;
Lp=Lpstar;
y=ystar;
n=nstar;
i=istar;
end;

steady;

check;

%% RANDOM SHOCKS BLOCK
shocks;
var epk; stderr 0.22;
var epv; stderr 0.21;
var er; stderr 0.17;
var ezf; stderr 0.03;
var ezp; stderr 0.04;
var epk, epv = 0.71;
var epk, er = 0.35;
var epk, ezf = -0.19;
var epk, ezp = -0.01;
var epv, er = 0.29;
var epv, ezf = -0.30;
var epv, ezp = -0.05;
var er, ezf = -0.05;
var er, ezp = -0.03;
var ezf, ezp = 0.06;
end;

stock_simul(order=1, hp_filter = 100,irf=50)y yf yp c i if ip n nf np nx v;




