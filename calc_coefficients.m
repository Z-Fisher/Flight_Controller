% MATAERO
%
% Copyright © 2011 by David K. Schmidt,
% published and distributed by McGraw-Hill with permission.
%
% This program computes the lift and moment coefficients for
% multiple wing surfaces. The tip boundary conditions of cl=0 is enforced.

NW=1; %Number of wings
N=[19]; %Number of segments subdividing wing
XLE=[-0.2793 -0.1969 -0.1723 -0.1477 -0.1231 -0.0985 -0.0739 -0.0492 -0.0246 0 -0.0246 -0.0492 -0.0739 -0.0985 -0.1231 -0.1477 -0.1723 -0.1969 -0.2793]; %Leading Edge x-coord for each segment (dim is value of N)
YLE=[-0.8268 -0.7366 -0.6431 -0.5512 -0.4594 -0.3675 -0.2756 -0.1837 -0.0919 0 0.0919 0.1837 0.2756 0.3675  0.4594 0.5512 0.6431 0.7366 0.8268]; %LE y-coord along wing
CHORD=[0 0.2924 0.4408 0.5467 0.6352 0.7104 0.7662 0.7990 0.6645 0.6599 0.6645 0.7990 0.7662 0.7104 0.6352 0.5467 0.4408 0.2924 0];
ih=[0]; %dim equal to val of NW
atwst=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
gtwst=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
refwng=1;
alpha=3;
AOA=[-6 -4 -2 0 2 4 6 8];
Cl=[-0.65 -0.45 -0.25 0 0.2 0.45 0.65 0.8];
Cm=[0 0 0 0 0 0 0 0];

clear G S MAC CL CM MIA
%
% Calculate the 1/4 chord location of each wing section
% and initialize the circulation distribution
%
k=0;
for i=1:NW,
  for j=1:N(i),
    k=k+1;
    xqc(k)=XLE(k)+0.25*CHORD(k);
    aeff=alpha+ih(i)+gtwst(k);
    sctn=atwst(k);
    [cl,cm]=gtclcm(sctn,aeff,AOA,Cl,Cm);
    gamma(k,1)=0.5*CHORD(k)*cl; % Does not include Vinfinity
  end
end
G=gamma;
clear k i j aeff sctn cl cm
%
% Get the Influence Coefficient Matrix
%
AIC=gtAIC(NW,N,xqc,YLE);
%
% Get the zero lift angles of attack for the tip sections
% and initialize the offset distances for the tip vorticies
%
k=0;
indx=0;
for i=1:NW,
  k=k+1;
  indx=indx+1;
  aol(indx)=fndzlaoa(atwst(k),Cl,AOA);
  offst(indx)=0.0;
  k=k+N(i)-1;
  indx=indx+1;
  aol(indx)=fndzlaoa(atwst(k),Cl,AOA);
  offst(indx)=0.0;
end
clear k indx i
%
% Begin the interation process
%
delta=.7;
cmax=7;
for count=2:cmax,
%
% set offset values such that cl=0 at the tips
  offst=tipbc(NW,N,xqc,YLE,gamma,offst,AIC,alpha,ih,gtwst,aol);
%
% get the induced angles of attack
  aind=gtinducd3(NW,N,xqc,YLE,gamma,offst,AIC);
%
% calculate the new circulation distribution
  k=0;
  for i=1:NW,
    for j=1:N(i),
      k=k+1;
      aeff=alpha+ih(i)+gtwst(k)+aind(k);
      sctn=atwst(k);
      [cl(k),cm(k)]=gtclcm(sctn,aeff,AOA,Cl,Cm);
      newgam=0.5*CHORD(k)*cl(k);
      gamma(k,1)=gamma(k,1)+delta*(newgam-gamma(k,1));
    end
  end
  G(:,count)=gamma;
%
% repeate this procedure
end
%cl
%cm
%aind
clear xqc AIC aol gamma count offst k i j aeff sctn newgam
%
% Now, calculate Lift and Moment coefficients
%
[S,MAC,XAC,CL,CM,MIA]=gtCLaCM(NW,N,XLE,YLE,CHORD,cl,cm,aind);
Sref=S(refwng);
MACref=MAC(refwng);
for i=1:NW,
	CMac(i)=CM(i)+CL(i)*XAC(i)/MAC(i);
end
CLtotal=0.0;
CMtotal=0.0;
CM0tot=0.0;
for i=1:NW,
  CLtotal=CLtotal+CL(i)*(S(i)/Sref);  
  CMtotal=CMtotal+CM(i)*(S(i)/Sref)*(MAC(i)/MACref);
  CM0tot=CM0tot+CMac(i)*(S(i)/Sref)*(MAC(i)/MACref);
end
clear i
S
Sref
MAC
MACref
XAC
%XACref
CL
CLtotal
CM
CMtotal
CMac
CM0tot
MIA