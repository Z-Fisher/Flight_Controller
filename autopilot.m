clc; clear all;
format long g

%Aircraft Specs
g=32.2; %accel of grav (ft/s^2)
b_w=1.635; %wing span 0.504 (m) or 1.6535 ft
sweep_LE=15; %Leading edge sweep (deg)
sweep_c2=-3.51; %Half chord sweep (deg)
sweep_vc2=7.79; %tail fin half chord sweep (deg)
alph_wing=3; %AOA of wing (deg)
alph_0wing=0; %AOA of zero lift coefficient
delX_AC=0.0709; % X coordinate of AC (ft)
Y_MAC=0.367*3.28; %Y coordinate of MAC (ft), 0.367 meters
A=2.856; %Aspect Ratio of wing
A_v=8.177; %Aspect ratio of vertical fin
c_lalph=0.1921*(180/pi); %section lift coefficient (1/rad)
a=343*3.28; %speed of sound (ft/s)
u=15.2*3.28; %Cruise speed (ft/s)
M=u/a; %Mach number
B=sqrt(1-M^2);
k=c_lalph/(2*pi);
U_0=49.9; %free stream velocity (ft/s)
q_inf=0.5*0.002378*U_0^2; %dynamic pressure (lb/ft^2)
q_h=0.9*q_inf; %dynamic pressure at tail (horizontal)
S_w=0.9534; %Wing surface area (ft^2)
S_v=0.0107; %Total vertical fin area (ft^2)
S_ail=0.1044; %total aileron surface area (ft^2)
m=0.558; %mass of vehicle (lbs)
I_xx=0.0007519; %slug*ft^2
I_yy=0.00117; %slug*ft^2
I_zz=0.00188; %slug*ft^2

%Effectiveness Calculations
e=0.75; %Oswald's span-efficient factor 'Assumption'
phi_T=0; %pg 257
alph_0=alph_0wing;
cbar_w=0.6501; %MAC (ft)
d_T=0; %pg 257
x_t=1; %DN
T_0=0;% don't know. pg 257

C_Lalph=2.413065;
C_Lroll=-.000001;
C_Lbeta=0;
C_Lalphdot=0;
C_Lq=0; %pg 286
C_LdelE=0.355; %Elevator lift effectiveness /rad
C_L0=0; 
C_Lu=0; 
C_Lp=0; 
C_Lr=0; 
C_LdelA=-0.001664; %vehicle's aileron rolling-moment effectiveness
C_LdelR=0; %vehicle's rudder rolling-moment effectiveness
C_L=0.8; 

C_D0=0.16; %Parasite drag coeff of vehicle
C_D=0.2452; %Vehicle drag coefficient
C_Du=0; %used in X_u equation
C_Dalph=2*C_Lalph^2/(pi*A*e); %Vehicle AOA drag effectiveness
C_Dalphdot=0;
C_Dq=0; %pg 286
C_DdelE=0; %vehicle elvator drag effectiveness

C_PXu=0; %pg 269
C_PX0=0; 
C_PZu=0; %pg 269
C_PZ0=0; 
C_PMu=0; %pg 270
C_PM0=0;
C_PMalph=0; 

C_Sbetav=-2*pi*A_v/(2+sqrt((A_v^2)*(B^2)*(1+(tan(sweep_vc2)^2)/B^2)+4));
C_Sbeta=C_Sbetav*q_h*S_v/(q_inf*S_w); %Vehicle sideslip side-force effectiveness
C_Sp=0; %Steady unaccelerated flight
C_Sr=0;
C_SdelA=0; %Vehicle aileron side-force effectiveness
C_SdelR=0; %Vehicle rudder side-force effectiveness

C_Mu=0; 
C_M0=0; %steady unaccelerated flight
C_Malph=0.006048; %vehicle AOA pitching-moment effectiveness
C_Malphdot=-0.0706;
C_Mq=-0.660064; %pg 287
C_MdelE=-0.87; %vehicle elevator pitching-moment effectiveness /rad

C_Nbeta=0.070381; %Vehicle sideslip yawing moment effectiveness
C_N0=0; %steady unaccelerated flight
C_Np=-0.049271;
C_Nr=-0.032935;
C_NdelA=0.000169; %pg 209. vehicle aileron yawing moment effectiveness
C_NdelR=0; 

%Force Equations
X_u__X_Pu=(q_inf*S_w/m)*(-1*(C_Du+2*C_D0/U_0)+(C_PXu+2*C_PX0/U_0));
X_alph=-(q_inf*S_w/m)*(-C_Dalph+C_L0);
X_alphdot=-q_inf*S_w*C_Dalphdot/m;
X_q=-q_inf*S_w*C_Dq/m;
X_delE=-q_inf*S_w*C_DdelE/m;
X_T=-cos(phi_T+alph_0)/m;

Y_beta=q_inf*S_w*C_Sbeta/m;
Y_p=q_inf*S_w*C_Sp/m;
Y_r=q_inf*S_w*C_Sr/m;
Y_delA=q_inf*S_w*C_SdelA/m;
Y_delR=q_inf*S_w*C_SdelR/m;

Z_u__Z_Pu=-0.3697;
Z_alph=(-q_inf*S_w/m)*(C_Lalph+C_D0);
Z_alphdot=-q_inf*S_w*C_Lalphdot/m; %zero
Z_q=-q_inf*S_w*C_Lq/m;
Z_delE=-q_inf*S_w*C_LdelE/m;
Z_T=-sin(phi_T+alph_0)/m;

%Moment Equations
L_beta=-11;%(q_inf*S_w*b_w/I_xx)*C_Lbeta;
L_p=(q_inf*S_w*b_w/I_xx)*C_Lp;
L_r=(q_inf*S_w*b_w/I_xx)*C_Lr;
L_delA=(q_inf*S_w*b_w/I_xx)*C_LdelA;
L_delR=(q_inf*S_w*b_w/I_xx)*C_LdelR;

M_u__M_Pu=(q_inf*S_w*cbar_w/I_yy)*(C_Mu+2*C_M0/U_0)+(C_PMu+2*C_PM0/U_0);
M_alph=(q_inf*S_w*cbar_w/I_yy)*C_Malph;
M_Palph=(q_inf*S_w*cbar_w/I_yy)*C_PMalph;
M_alphdot=(q_inf*S_w*cbar_w/I_yy)*C_Malphdot; %nonzero
M_q=(q_inf*S_w*cbar_w/I_yy)*C_Mq;
M_delE=(q_inf*S_w*cbar_w/I_yy)*C_MdelE;
M_T=(d_T*cos(phi_T)-x_t*sin(phi_T))/I_yy;

N_beta=(q_inf*S_w*b_w/I_zz)*C_Nbeta;
N_p=(q_inf*S_w*b_w/I_zz)*C_Np;
N_r=(q_inf*S_w*b_w/I_zz)*C_Nr;
N_delA=(q_inf*S_w*b_w/I_zz)*C_NdelA;
N_delR=(q_inf*S_w*b_w/I_zz)*C_NdelR;

%Longitudinal Matrices Calculations
A_long=[X_u__X_Pu+(X_alphdot*Z_u__Z_Pu/(U_0-Z_alphdot)) X_alph+(X_alphdot*Z_alph/(U_0-Z_alphdot)) -g X_q+X_alphdot*(U_0+Z_q)/(U_0-Z_alphdot) 0;
    (Z_u__Z_Pu)/(U_0-Z_alphdot) Z_alph/(U_0-Z_alphdot) 0 (U_0+Z_q)/(U_0-Z_alphdot) 0;
    0 0 0 1 0;
    M_u__M_Pu+(M_alphdot*(Z_u__Z_Pu)/(U_0-Z_alphdot)) M_alph+M_Palph+M_alphdot*Z_alph/(U_0-Z_alphdot) 0 M_q+M_alphdot*(U_0+Z_q)/(U_0-Z_alphdot) 0;
    0 -U_0 U_0 0 0]
B_long=[X_delE+X_alphdot*Z_delE/(U_0-Z_alphdot) X_T+X_alphdot*Z_T/(U_0-Z_alphdot);
    Z_delE/(U_0-Z_alphdot) Z_T/(U_0-Z_alphdot);
    0 0;
    M_delE+M_alphdot*Z_delE/(U_0-Z_alphdot) M_T+M_alphdot*Z_T/(U_0-Z_alphdot);
    0 0]
C=eye(5);
D=zeros(5,2);

%Lateral Matrices Calculations
A_lat=[Y_beta/U_0 g/U_0 Y_p/U_0 (Y_r/U_0)-1 0;
    0 0 1 0 0;
    L_beta 0 L_p L_r 0;
    N_beta 0 N_p N_r 0;
    0 0 0 1 0]
B_lat=[Y_delA/U_0 Y_delR/U_0;
    0 0;
    L_delA L_delR;
    N_delA N_delR;
    0 0]

%S.S. Systems
sys_lat=ss(A_lat,B_lat,C,D); %LATERAL
sys_long=ss(A_long,B_long,C,D); %LONGITUDINAL

%Dampers
Kq=-2; %longitudinal
A_long_aug=A_long-B_long(:,1)*Kq*C(4,:);
sys_long_aug=ss(A_long_aug,B_long,C,D);

Kr=0.18; %lateral
A_lat_aug=A_lat-B_lat(:,1)*Kr*C(4,:);
sys_lat_aug=ss(A_lat_aug,B_lat,C,D);

%Aileron input
t=[0:.1:10]; %time interval
u_aileron=zeros(2,101);
u_aileron(1,1:10)=-1*(pi/180); %1 sec by 1 degree

%Elevator Doublet
u_elevator=zeros(2,101);
u_elevator(1,1:10)=-1*(pi/180); %1 sec at -1 degree
u_elevator(1,11:20)=1*(pi/180); %1 sec at +1 degree

%lsim plotting
figure(1);
lsim(sys_lat,u_aileron,t) %1 sec aileron input
figure(2);
lsim(sys_long,u_elevator,t) %2 sec elevator doublet
figure(3);
lsim(sys_lat_aug,u_aileron,t) %augmented lat
figure(4);
lsim(sys_long_aug,u_elevator,t) %augmented long

%Eigenvalues
eigs_lat=eig(A_lat)
eigs_long=eig(A_long)
eigs_lat_aug=eig(A_lat_aug)
eigs_long_aug=eig(A_long_aug)

%Autopilot Altitude Hold (longitudinal. pitch-attitude)
zsysauglong=zpk(sys_long_aug);
uelevaug_long=zsysauglong(1,1);
alphaelevaug_long=zsysauglong(2,1);
thetaelevaug_long=zsysauglong(3,1);
qelevaug_long=zsysauglong(4,1);
[z1,p1,k1]=zpkdata(alphaelevaug_long); 
[z2,p2,k2]=zpkdata(thetaelevaug_long);
gammatheta=1-zpk(z1,z2,k1/k2);
hgamma=tf(U_0,[1 0]);
sysgam=series(gammatheta,hgamma);
PI=tf([1 3],[1 0]); %PI Compensator
LC=tf([1 0.0001],[1 1]); %Lead Compensator
Kh=-0.00019; %Gain
sysfinal=Kh*PI*sysgam*LC; %Form Open Loop System for Altitude
figure(5);
margin(sysfinal);
[Afin,Bfin,Cfin,Dfin]=ssdata(sysfinal);
Afincl=Afin-Bfin*1*Cfin;
sysclalt=ss(Afincl,Bfin,Cfin,Dfin);
zsysclalt=zpk(sysclalt);

%Autopilot Heading Hold (lateral. bank-angle)
syscl=sys_lat_aug;
hdg=tf(g/U_0,[1 0]);
syspsi=series(syscl,hdg,[2],[1]);
[Apsi,Bpsi,Cpsi,Dpsi]=ssdata(syspsi);
Kpsi=-0.05;
Bpsi=Bpsi(:,1);
Aclpsi=Apsi-Bpsi*Kpsi*Cpsi;
Bclpsi=Bpsi*Kpsi;
Cclpsi=C;
Cclpsi=[[0;0;0;0;0] C];
Cclpsi(6,:)=Cpsi;
Dclpsi=[0;0;0;0;0;0]; 
sysclpsi=ss(Aclpsi,Bclpsi,Cclpsi,Dclpsi);
zsysclpsi=zpk(sysclpsi);
hdgcmd=zsysclpsi(6,1);

%Step Change of 100ft in altitude and 3deg in heading
figure(7);
step(100*sysclalt,100) %altitude change response
figure(8);
step(0.05236*sysclpsi,100) %heading response (3deg in radians)
figure(9);
step(0.05236*hdgcmd,100)
