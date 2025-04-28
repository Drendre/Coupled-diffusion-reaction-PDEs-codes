%% Diffusion with simple coupling code
clear;
%% Parameters
S=15; %number of different time step sizes
D=1; A=-0.6; alpha=0.87; gamma=A+alpha+1; %% coefficients
c1=1; c2=0; c3=0; c4=0; a=4; b=9; %% coefficients
ti=0.04; TIME=0.1; tf=ti+TIME;  %% initial time, Length of time interval, final time
x0=-2; Lx=4; xf=Lx+x0; %Space domain
Nx=200; Nx=Nx+1; Dx=(xf-x0)/(Nx-1); %Space discretization, Nx must be odd
xaxis=zeros(Nx,1); %space axis
for i=1:Nx 
    xaxis(i)=x0+(i-1)*Dx; 
end 
u0=zeros(Nx,1);uex=zeros(Nx,1);w0=zeros(Nx,1);wex=zeros(Nx,1); %%Initial and exact solution
MaxD=zeros(S,10);L2D=zeros(S,10);Axhstep=zeros(S,1); %errors and stepsize-axis
p=1-gamma/2-alpha/2;ee=zeros(Nx,1);ee1=zeros(Nx,1);ee2=zeros(Nx,1);ee13=zeros(Nx,1);ee23=zeros(Nx,1);e1e=zeros(Nx,1);e1e2=zeros(Nx,1);
AbsDU_UQ= zeros(Nx,1,'double');AbsDU_UQ1= zeros(Nx,1,'double');AbsDU_ULH= zeros(Nx,1,'double');AbsDU_ULH1= zeros(Nx,1,'double');AbsDU_UADE1= zeros(Nx,1,'double');AbsDU_UADE= zeros(Nx,1,'double');AbsDU_Uuodd1= zeros(Nx,1,'double'); AbsDU_UDF= zeros(Nx,1,'double');
AbsDU_UCCL= zeros(Nx,1,'double');AbsDU_UCCL1= zeros(Nx,1,'double');AbsDU_UCpC1= zeros(Nx,1,'double');AbsDU_UCpC= zeros(Nx,1,'double');AbsDU_UPI1= zeros(Nx,1,'double');AbsDU_UPI= zeros(Nx,1,'double');AbsDU_Uuodd= zeros(Nx,1,'double');AbsDU_UDF1= zeros(Nx,1,'double');
AveD=zeros(S,1,'double');ARE=zeros(S,1,'double');Sum_Log_MaxD=zeros(S,1,'double'); Sum_Log_AveD=zeros(S,1,'double');Sum_Log_DifQ=zeros(S,1,'double');  
q=(sqrt(4*a*b-2*alpha*gamma+gamma^2+alpha^2))/2;
f1=gamma-alpha+2*q; % f2=gamma-alpha-2*q;
%%%%%%%%%Initial condition
for i=1:Nx
kk=xaxis(i)^2/ti; 
u0(i)=xaxis(i)/ti^(1/2)/(ti^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)); %+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
w0(i)=xaxis(i)/ti^(1/2)/(ti^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
%%%Final analitical solution
kk=xaxis(i)^2/tf;
uex(i)=xaxis(i)/tf^(1/2)/(tf^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p-q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p-q,3/2,kk/4));
wex(i)=xaxis(i)/tf^(1/2)/(tf^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
end
NN=[max(abs(u0)) max(abs(w0)) max(abs(uex)) max(abs(wex))]; % Maximum norm
Inorm=1/max(NN); % Inverse maximum norm
u0=Inorm*u0;w0=Inorm*w0;uex=Inorm*uex;wex=Inorm*wex; % Normalization
% plot(xaxis,u0,xaxis,uex,xaxis,w0,':s',xaxis,wex,':x')
%%%%%%%%%%%%%%
w=zeros(Nx,1);w=zeros(Nx,1);u=zeros(Nx,1);a0=zeros(Nx,1);a01=zeros(Nx,1); ULH=zeros(Nx,1);ULH1=zeros(Nx,1); UCpC=zeros(Nx,1);UCpC1=zeros(Nx,1);
uTEMP=zeros(Nx,1);wTEMP=zeros(Nx,1);uTEMP1=zeros(Nx,1);wTEMP1=zeros(Nx,1);uTEMP=zeros(Nx,1);wTEMP=zeros(Nx,1);uTEMP2=zeros(Nx,1);wTEMP2=zeros(Nx,1);%%Values of U, U2 is the latter
E=D/Dx^2; r2=1/Dx/2; %Matrix elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ih=1:S  %% Big loop for the timestep, ih is the Serial Number
ih
h=TIME/2^1/(2^ih);  T = ceil(TIME/h-0.1)+1; %(new, decreased) time step size and (increased) number of required time-steps
taxis=zeros(T+1,1); %%physical time
for t=1:T+1
      taxis(t)=ti+(t-1)*h;
end
Axhstep(ih)=h; %%Timestep axis
r=h*E;
eed=exp(-2*r);emid=exp(-r);e13=exp(-2*r/3);e23=exp(-4*r/3);
  if(r>0.01)  
      e1ed=1-eed;e1mid=1-emid;e13e=1-e13;e23e=1-e23;
  else
    e1mid=0; e1ed=0; e13e=0;e23e=0;
    for s=1:10
        sf=factorial(s);
      e1mid=e1mid-(-r)^s/sf; e1ed=e1ed-(-2*r)^s/sf;
    e13e=e13e-(-2*r/3)^s/sf; e23e=e23e-(-2*r*2/3)^s/sf;
    end
  end
F=1-e1ed/2/r;Fmid=(1-e1mid/(r));
 if(h>0.00001 && h<0.0001) rfun(ih)=2000000000*h*h; end %
%%%Boundary conditions left and right
bound1u=zeros(T,1);bound1w=zeros(T,1);boundNu=zeros(T,1);boundNw=zeros(T,1);
boundNu12=zeros(T,1);boundNw12=zeros(T,1);boundNu13=zeros(T,1);boundNw13=zeros(T,1);boundNu23=zeros(T,1);boundNw23=zeros(T,1);
bound1u12=zeros(T,1);bound1w12=zeros(T,1);bound1u13=zeros(T,1);bound1w13=zeros(T,1);bound1u23=zeros(T,1);bound1w23=zeros(T,1);
tic
for tind=1:T %%tind is the index of the time, always integer
    t=ti+tind*h;  %t is the physical time in the system
    kk=x0^2/t;
bound1u(tind)=Inorm*x0/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
bound1w(tind)=Inorm*x0/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
   kk=xf^2/t; 
boundNu(tind)=Inorm*xf/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
boundNw(tind)=Inorm*xf/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
t=ti+(tind-1/2)*h;  kk=x0^2/t; 
bound1u12(tind)=Inorm*x0/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
bound1w12(tind)=Inorm*x0/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
kk=xf^2/t; 
boundNu12(tind)=Inorm*xf/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
boundNw12(tind)=Inorm*xf/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
t=ti+(tind-2/3)*h;  kk=x0^2/t; 
bound1u13(tind)=Inorm*x0/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
bound1w13(tind)=Inorm*x0/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
kk=xf^2/t; 
boundNu13(tind)=Inorm*xf/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
boundNw13(tind)=Inorm*xf/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
t=ti+(tind-1/3)*h;  kk=x0^2/t; 
bound1u23(tind)=Inorm*x0/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
bound1w23(tind)=Inorm*x0/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
kk=xf^2/t; 
boundNu23(tind)=Inorm*xf/t^(1/2)/(t^alpha)*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4));%+c2*hypergeom(p+q,3/2,kk/4)+c3*kummerU(p-q,3/2,kk/4)+c4*kummerU(p+q,3/2,kk/4));
boundNw23(tind)=Inorm*xf/t^(1/2)/(t^gamma)/2/a*exp(-kk/4)*(c1*hypergeom(p-q,3/2,kk/4)*f1);%+c2*hypergeom(p+q,3/2,kk/4)*f2+c3*kummerU(p-q,3/2,kk/4)*f1+c4*kummerU(p+q,3/2,kk/4)*f2);
end
toc

%% The numerical methods %%%%%%%%%%%%%%%% 
%% FTCS method 
u=u0; w=w0;
%% Big loop for time
for tind=1:T-1 %%tind is the index of the time, always integer
    t=ti+(tind-1/2)*h;
  for i=2:Nx-1
        uTEMP(i)=u(i)+r*(u(i-1)+u(i+1)-2*u(i)) + h*a*w(i)*t^A; %%% u variable
        wTEMP(i)=w(i)+r*(w(i-1)+w(i+1)-2*w(i)) + h*b*u(i)*t^(-A-2); %%% w variable
  end
    uTEMP(1)=bound1u(tind); uTEMP(Nx)=boundNu(tind); u=uTEMP; 
    wTEMP(1)=bound1w(tind); wTEMP(Nx)=boundNw(tind); w=wTEMP;
end %% End big loop for time

err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,1)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,2)=err_w;

Diff=abs(uex-u); % maximum error
if(Diff<1000) MaxD(ih,1)=max(Diff); end
Diff=abs(wex-w);
if(Diff<1000) MaxD(ih,2)=max(Diff); end
UEE=u;

%% Adam-Bashforth method, AB2
 u=u0;w=w0; % these will be the latest values
 u2=u;w2=w; % these will be the old values
 tic
 t=ti+1/2*h;
for i=2:Nx-1 %%% zero-th stage, the latest/new values are updated
u(i)=u(i)*(1-2*r)+(u(i-1)+u(i+1))*r+h*a*w(i)*t^A; 
 w(i)=w(i)*(1-2*r)+(w(i-1)+w(i+1))*r+h*b*u(i)*t^(-A-2);
end
 u(1)=bound1u(1); u(Nx)=boundNu(1);  
 w(1)=bound1w(1); w(Nx)=boundNw(1);
 
 %% Big loop for time
for tind=2:T-1 %%tind is the index of the time, always integer
    t=ti+(tind-1/2)*h;
f10=zeros(Nx,1);f20=zeros(Nx,1); f11=zeros(Nx,1);f21=zeros(Nx,1);
for i=2:Nx-1   %%% first stage 
   f10(i)=u(i)*(-2*r)+(u(i-1)+u(i+1))*r +h*a*w(i)*t^A; %% using the new values
   f11(i)=w(i)*(-2*r)+(w(i-1)+w(i+1))*r +h*b*u(i)*t^(-A-2); %% using the new values
   f20(i)=u2(i)*(-2*r)+(u2(i-1)+u2(i+1))*r +h*a*w2(i)*t^A; %% using the old values 
   f21(i)=w2(i)*(-2*r)+(w2(i-1)+w2(i+1))*r +h*b*u2(i)*t^(-A-2); %% using the old values
 end
 uTEMP=u+(3*f10-f20)/2; wTEMP=w+(3*f11-f21)/2;
  uTEMP(1)=bound1u(tind); uTEMP(Nx)=boundNu(tind);  
    wTEMP(1)=bound1w(tind); wTEMP(Nx)=boundNw(tind);
 u2=u; u=uTEMP;w2=w;w=wTEMP;
end
toc
UEE=u;
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,3)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,4)=err_w;

Diff=abs(uex-u); %%% maximum error
if(Diff<1000) MaxD(ih,3)=max(Diff); end
Diff=abs(wex-w);
if(Diff<1000) MaxD(ih,4)=max(Diff); end
%% ADE method 
u=u0;w=w0;tic
% t=ti+1/2*h;
for tind=1:T-1 %%tind is the index of the time, always integer
    t=ti+(tind-1/2)*h; uTEMP=u;uTEMP1=u;wTEMP=w;wTEMP1=w;
 uTEMP(1)=bound1u(tind); wTEMP(1)=bound1w(tind);
for i=2:Nx-1 %%% from left to right
 a0(i)=(uTEMP(i-1)+uTEMP(i+1)-uTEMP(i))*r+(h*a*w(i)*t^A);
 a01(i)=(wTEMP(i-1)+wTEMP(i+1)-wTEMP(i))*r+(h*b*u(i)*t^(-A-2));
 uTEMP(i)=(uTEMP(i)+a0(i))/(1+r);
 wTEMP(i)=(wTEMP(i)+a01(i))/(1+r);
end
 uTEMP1(Nx)=boundNu(tind); wTEMP1(Nx)=boundNw(tind);
for i=2:Nx-1   %%% from right to left  
  j=Nx-i+1 ;
  au(j)=(uTEMP1(j-1)+uTEMP1(j+1)-uTEMP1(j))*r+h*a*w(j)*t^A;
  aw(j)=(wTEMP1(j-1)+wTEMP1(j+1)-wTEMP1(j))*r+h*b*u(j)*t^(-A-2);
 uTEMP1(j)=(uTEMP1(j)+au(j))/(1+r);
 wTEMP1(j)=(wTEMP1(j)+aw(j))/(1+r);
end
 uTEMP(Nx)=boundNu(tind);wTEMP(Nx)=boundNw(tind);
 uTEMP1(1)=bound1u(tind);wTEMP1(1)=bound1w(tind);
u=(uTEMP+uTEMP1)/2;w=(wTEMP+wTEMP1)/2;
end 
toc
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,5)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,6)=err_w;


AbsDU_UADE=max(abs(u-uex));MaxD(ih,5)=max(AbsDU_UADE);%%%error
AbsDU_UADE1=max(abs(w-wex));MaxD(ih,6)=max(AbsDU_UADE1);%%%error

%%  Dufort-Frankel method, DF
 u=u0;uold=u;w=w0;wold=w;  t=ti+h/2;
tic
  for i=2:Nx-1  %%% zero-th stage, UPFD formula
       u(i)=(u(i)+r*(u(i-1)+u(i+1))+h*a*w(i)*t^A)/(1+2*r); 
      w(i)=(w(i)+r*(w(i-1)+w(i+1))+h*b*uold(i)*t^(-A-2))/(1+2*r);
  end
   u(1)=bound1u(1); u(Nx)=boundNu(1);  
    w(1)=bound1w(1); w(Nx)=boundNw(1);
 %% Big loop for time
for tind=2:T-1
    for i=2:Nx-1  %%% second stage 
        t=ti+(tind-1/2)*h;
         uTEMP(i)=((1-2*r)*uold(i)+2*r*(u(i-1)+u(i+1))+2*h*a*w(i)*t^A)/(1+2*r);
         wTEMP(i)=((1-2*r)*wold(i)+2*r*(w(i-1)+w(i+1))+2*h*b*u(i)*t^(-A-2))/(1+2*r);
     end
   uTEMP(1)=bound1u(tind); uTEMP(Nx)=boundNu(tind);  
   wTEMP(1)=bound1w(tind); wTEMP(Nx)=boundNw(tind);
   uold=u;wold=w;
   u=uTEMP;w=wTEMP;
end
toc
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,7)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,8)=err_w;

AbsDU_UDF=max(abs(u-uex));MaxD(ih,7)=max(AbsDU_UDF);%%% Max error
AbsDU_UDF1=max(abs(w-wex));MaxD(ih,8)=max(AbsDU_UDF1);%%% Max error

%% Original odd even hopscotch, OOEH
u=u0;w=w0;tic 
for i=1:1:Nx %%%
cpar(i)=(-1)^i;
end
for tind=1:T-1 %%%tind is the index of the time, always integer
    t=ti+(tind-1/2)*h; 
 for i=2:Nx-1   %%%First stage, Explicit euler 
     if (cpar(i)<0)  
     u(i)=u(i)*(1-2*r)+(u(i-1)+u(i+1))*r+h*a*w(i)*t^A; 
     w(i)=w(i)*(1-2*r)+(w(i-1)+w(i+1))*r+h*b*u(i)*t^(-A-2);
     end
 end
    u(1)=bound1u(tind); u(Nx)=boundNu(tind);  
    w(1)=bound1w(tind); w(Nx)=boundNw(tind);
for i=2:Nx-1 %%% second stage, BTCS formula
if (cpar(i)>0)
      u(i)=(u(i)+r*(u(i-1)+u(i+1))+h*a*w(i)*t^A)/(1+2*r); 
      w(i)=(w(i)+r*(w(i-1)+w(i+1))+h*b*u(i)*t^(-A-2))/(1+2*r);  
end
end
 cpar=-cpar;
end %%%End big loop for time
toc
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,9)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,10)=err_w;

AbsDU_Uuodd=max(abs(u-uex));MaxD(ih,9)=max(AbsDU_Uuodd);%%% Max error
AbsDU_Uuodd1=max(abs(w-wex));MaxD(ih,10)=max(AbsDU_Uuodd1);%%%error

%% CpC method
u=u0;w=w0; tic
for tind=1:T-1 %%tind is the index of the time, always integer
    t=ti+(tind-3/4)*h; 
for i=2:Nx-1   %%% first stage, predictor with half time step  
    a0=(u(i-1)+u(i+1))/2 +(h*a*w(i)*t^A)/(2*r);
    a01=(w(i-1)+w(i+1))/2 +(h*b*u(i)*t^(-A-2))/(2*r);
    uTEMP(i)=u(i)*emid+a0*e1mid; 
     wTEMP(i)=w(i)*emid+a01*e1mid;
end
 uTEMP(1)=bound1u12(tind); uTEMP(Nx)=boundNu12(tind);  
    wTEMP(1)=bound1w12(tind); wTEMP(Nx)=boundNw12(tind);
t=ti+(tind-1/2)*h;
for i=2:Nx-1 %%% second stage, corrector with full time step 
a0(i)=(uTEMP(i+1)+uTEMP(i-1))/2+(h*a*wTEMP(i)*t^A)/(2*r);
a01(i)=(wTEMP(i+1)+wTEMP(i-1))/2+(h*b*uTEMP(i)*t^(-A-2))/(2*r);
  u(i)=u(i)*eed+a0(i)*e1ed; 
  w(i)=w(i)*eed+a01(i)*e1ed;
end
u(1)=bound1u(tind); u(Nx)=boundNu(tind);  
 w(1)=bound1w(tind); w(Nx)=boundNw(tind);
end
toc
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,11)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,12)=err_w;

AbsDU_UCpC=max(abs(u-uex));MaxD(ih,11)=max(AbsDU_UCpC);%%%error
AbsDU_UCpC1=max(abs(w-wex));MaxD(ih,12)=max(AbsDU_UCpC1);%%%error

%% CCL method 
u=u0;w=w0; tic
for tind=1:T-1 %%tind is the index of the time, always integer
t=ti+(tind-5/6)*h;
for i=2:Nx-1 %%% first satge 1/3 length
a0(i)=(u(i-1)+u(i+1))/2+(h*a*w(i)*t^A)/(r*2);
a01(i)=(w(i-1)+w(i+1))/2+(h*b*u(i)*t^(-A-2))/(r*2);
uTEMP(i)=u(i)*e13+a0(i)*e13e; 
wTEMP(i)=w(i)*e13+a01(i)*e13e;
end
uTEMP(1)=bound1u13(tind); uTEMP(Nx)=boundNu13(tind); 
wTEMP(1)=bound1w13(tind); wTEMP(Nx)=boundNw13(tind);
t=ti+(tind-2/3)*h;
for i=2:Nx-1 %%% second stage 2/3 length
q1=(uTEMP(i-1)+uTEMP(i+1))/2+(h*a*wTEMP(i)*t^(A))/(r*2);
q2=(wTEMP(i-1)+wTEMP(i+1))/2+(h*b*uTEMP(i)*t^(-A-2))/(r*2);
uTEMP1(i)=u(i)*e23+q1*e23e; 
wTEMP1(i)=w(i)*e23+q2*e23e;
end
uTEMP1(1)=bound1u23(tind); uTEMP1(Nx)=boundNu23(tind); 
wTEMP1(1)=bound1w23(tind); wTEMP1(Nx)=boundNw23(tind);
t=ti+(tind-1/2)*h;
for i=2:Nx-1 %%% third stage, full lenght
q1=(uTEMP1(i-1)+uTEMP1(i+1))/2+(h*a*wTEMP1(i)*t^(A))/(r*2);
q2=(wTEMP1(i-1)+wTEMP1(i+1))/2+(h*b*uTEMP1(i)*t^(-A-2))/(r*2);
Su=3/2*(q1-a0(i));Sw=3/2*(q2-a01(i)); 
uTEMP(i)=u(i)*eed+a0(i)*e1ed+Su*F; 
wTEMP(i)=w(i)*eed+a01(i)*e1ed+Sw*F;
end
u=uTEMP;w=wTEMP;
u(1)=bound1u(tind); u(Nx)=boundNu(tind); 
w(1)=bound1w(tind); w(Nx)=boundNw(tind);
end 
toc
UC=u; WC=w;
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,13)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,14)=err_w;
AbsDU_UCCL=max(abs(u-uex));MaxD(ih,13)=max(AbsDU_UCCL);%%%error
AbsDU_UCCL1=max(abs(w-wex));MaxD(ih,14)=max(AbsDU_UCCL1);%%%error

%% Quadratic CLQQ method 
u=u0; u2=zeros(Nx,1); umid=zeros(Nx,1);u2mid=zeros(Nx,1);au0=zeros(Nx,1);
w=w0; w2=zeros(Nx,1); wmid=zeros(Nx,1);w2mid=zeros(Nx,1);aw0=zeros(Nx,1);
tic
for tind=1:T-1 %%tind is the index of the time, always integer
     t=ti+(tind-1/2)*h; tq=ti+(tind-3/4)*h;
   for i=2:Nx-1 %% first stage, CNe 
         au0(i)=(u(i-1)+u(i+1))/2 +(h*a*w(i)*t^(A))/(2*r);
         aw0(i)=(w(i-1)+w(i+1))/2 +(h*b*u(i)*t^(-A-2))/(2*r);
         uTEMP(i)=u(i)*eed+ au0(i)*e1ed;
         wTEMP(i)=w(i)*eed+ aw0(i)*e1ed;
   end
    uTEMP(1)=bound1u(tind); uTEMP(Nx)=boundNu(tind);
    wTEMP(1)=bound1w(tind); wTEMP(Nx)=boundNw(tind);
   for i=2:Nx-1  %%%Second stage, LNe 
      aunew=(uTEMP(i-1)+uTEMP(i+1))/2 +(h*a*wTEMP(i)*t^A)/(2*r);
      awnew=(wTEMP(i-1)+wTEMP(i+1))/2 +(h*b*uTEMP(i)*t^(-A-2))/(2*r);
      su=aunew-au0(i);sw=awnew-aw0(i);
   uTEMP2(i)= u(i)*eed + au0(i)*e1ed + su*F;
   wTEMP2(i)= w(i)*eed + aw0(i)*e1ed + sw*F;
   umid(i)= u(i)*emid + au0(i)*e1mid + su*Fmid/2;
   wmid(i)= w(i)*emid + aw0(i)*e1mid + sw*Fmid/2;   
   end
   umid(1)=bound1u12(tind); umid(Nx)=boundNu12(tind); %Boundary
   wmid(1)=bound1w12(tind); wmid(Nx)=boundNw12(tind);
   uTEMP2(1)=uTEMP(1);uTEMP2(Nx)=uTEMP(Nx);
   wTEMP2(1)=wTEMP(1);wTEMP2(Nx)=wTEMP(Nx);
    for ii=1:2 %%% automatic repetition of the third stage as the fourth stage
     for i=2:Nx-1   %%% third stage, Quadratic neighbour
         aunew=(uTEMP2(i-1)+uTEMP2(i+1))/2 +(h*a*wTEMP2(i)*t^A)/(2*r);
         awnew=(wTEMP2(i-1)+wTEMP2(i+1))/2 +(h*b*uTEMP2(i)*t^(-A-2))/(2*r);
         aumid=(umid(i-1)+umid(i+1))/2 +(h*a*wTEMP2(i)*tq^A)/(2*r);
         awmid=(wmid(i-1)+wmid(i+1))/2 +(h*b*uTEMP2(i)*tq^(-A-2))/(2*r);
         au00=au0(i); su=4*aumid-aunew-3*au00; gu=2*(aunew-2*aumid+au00);
         aw00=aw0(i); sw=4*awmid-awnew-3*aw00; gw=2*(awnew-2*awmid+aw00);
       u2(i)=u(i)*eed+ e1ed*(gu/2/r^2-su/2/r+au00)+ gu -gu/r +su ;
       w2(i)=w(i)*eed+ e1ed*(gw/2/r^2-sw/2/r+aw00)+ gw -gw/r +sw ;
        u2mid(i)=u(i)*emid+ e1mid*(gu/2/r^2-su/2/r+au00)+ gu/4-gu/2/r +su/2;
        w2mid(i)=w(i)*emid+ e1mid*(gw/2/r^2-sw/2/r+aw00)+ gw/4-gw/2/r +sw/2;
     end
     u2mid(1)=umid(1); u2mid(Nx)=umid(Nx); 
     w2mid(1)=wmid(1); w2mid(Nx)=wmid(Nx);
    u2(1)=uTEMP(1);u2(Nx)=uTEMP(Nx);
    w2(1)=wTEMP(1);w2(Nx)=wTEMP(Nx); %Boundary
    uTEMP2=u2;umid=u2mid; wTEMP2=w2;wmid=w2mid;
   end
  u=uTEMP2;w=wTEMP2;
end
toc
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,15)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,16)=err_w;
 AbsDU_UQ=max(abs(u-uex));MaxD(ih,15)=max(AbsDU_UQ);%%%error
AbsDU_UQ1=max(abs(w-wex));MaxD(ih,16)=max(AbsDU_UQ1);%%%error

%% Pseudo-implicit method, PI
u=u0;w=w0;tic
% t=ti+1/2*h;
for tind=1:T-1 %%tind is the index of the time, always integer
    t=ti+(tind-1/2)*h;
for i=2:Nx-1 %%% first stage, half time step 
a0(i)=(u(i-1)+u(i+1))*r/2+(h/2*a*w(i)*t^A);
a01(i)=(w(i-1)+w(i+1))*r/2+(h/2*b*u(i)*t^(-A-2));
 uTEMP(i)=(u(i)+a0(i))/(1+r); 
 wTEMP(i)=(w(i)+a01(i))/(1+r);
end
uTEMP(1)=bound1u12(tind); uTEMP(Nx)=boundNu12(tind);  
 wTEMP(1)=bound1w12(tind); wTEMP(Nx)=boundNw12(tind);
for i=2:Nx-1   %%% second stage, full time step  
  au(i)=(uTEMP(i-1)+uTEMP(i+1))*r+h*a*wTEMP(i)*t^A;
  aw(i)=(wTEMP(i-1)+wTEMP(i+1))*r+h*b*uTEMP(i)*t^(-A-2);
  u(i)=(u(i)*(1-r)+au(i))/(1+r);  
  w(i)=(w(i)*(1-r)+aw(i))/(1+r);
 end
u(1)=bound1u(tind); u(Nx)=boundNu(tind);  
w(1)=bound1w(tind); w(Nx)=boundNw(tind);
 end 
toc
upi=u; 
err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,17)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,18)=err_w;
AbsDU_UPI=max(abs(u-uex));MaxD(ih,17)=max(AbsDU_UPI);%%% maximum error
AbsDU_UPI1=max(abs(w-wex));MaxD(ih,18)=max(AbsDU_UPI1);

%% Leapfrog-Hopscotch method, LH
u=u0; w=w0;  
tic
t=ti+h/4;
 for i=2:2:Nx-1  %%% zero-th stage for even nodes
      u(i)=(u(i)+r/2*(u(i-1)+u(i+1))+h/2*a*w(i)*t^A)/(1+r); 
      w(i)=(w(i)+r/2*(w(i-1)+w(i+1))+h/2*b*u(i)*(t^(-A-2)))/(1+r);
 end
%% Big loop for time 
for tind=1:T-1 %%tind is the index of the time, always integer
   t=ti+(tind-1/2)*h;  %t is the physical time in the system
    u(1)=bound1u(tind); u(Nx)=boundNu(tind);  
    w(1)=bound1w(tind); w(Nx)=boundNw(tind);
    uTEMP=u;
  for i=3:2:Nx-1  %%%  odd stage  
     u(i)=(r*(u(i-1)+u(i+1))+(1-r)*u(i)+h*a*w(i)*t^A)/(1+r);
     w(i)=(r*(w(i-1)+w(i+1))+(1-r)*w(i)+h*b*uTEMP(i)*t^(-A-2))/(1+r);
  end
  t=t+h/2;uTEMP=u;
  for i=2:2:Nx-1  %%% even stage
    if (taxis(tind+1)<tf-0.8*h) %% if we are far from the end, full time step
        u(i)=(r*(u(i-1)+u(i+1))+(1-r)*u(i)+h*a*w(i)*t^A)/(1+r);
        w(i)=(r*(w(i-1)+w(i+1))+(1-r)*w(i)+h*b*uTEMP(i)*t^(-A-2))/(1+r);
    else %% if we are close to the end, half time step 
        u(i)=(r/2*(u(i-1)+u(i+1))+(1-r/2)*u(i)+h/2*a*w(i)*t^A)/(1+r/2);
        w(i)=(r/2*(w(i-1)+w(i+1))+(1-r/2)*w(i)+h/2*b*uTEMP(i)*t^(-A-2))/(1+r/2);
            end
  end
end
toc
ulh=u;

err_u=sqrt(Dx*sum((u-uex).^2)); %L2 error
L2D(ih,19)=err_u;
err_w=sqrt(Dx*sum((w-wex).^2)); %L2 error
L2D(ih,20)=err_w;

AbsDU_ULH=max(abs(u-uex));% errors
 MaxD(ih,19)=max(AbsDU_ULH); 
AbsDU_ULH1=max(abs(w-wex));
MaxD(ih,20)=max(AbsDU_ULH1); 


end %%% end of big loop for the set of time step sizes

%% Aggregated errors
for ie=1:20   
Sum_Log_MaxD(ie)=double(sum(log10(MaxD(:,ie))));
ARE(ie)=double((Sum_Log_MaxD(ie)));%+Sum_Log_AveD(ie)+Sum_Log_DifQ(ie))/(3.0 )); 
end

%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%

Blue=[0 , 0, 0.95];Red	= [1, 0, 0.1];
greenblue	 =    [0.4,	0.8,	0.7];orange = [1 0.5 0];
purple = [0.5 , 0, 0.5];purple2 = [0.6 , 0, 0.4];
darkgreen = [0.08, 0.7, 0.1];darkgreen2 = [0.02, 0.75, 0.08];
chocolatedark= [0.068, 0.79, 0.45];
Black	= [0, 0, 0];
magenta= [1, 0.4, 1]; darkredmagenta  = [0.9, 0.3, 1];
yellow  = [0.9 0.9 0];yellowdark = [0.95 0.8 0.01];
black = [0, 0, 0];  black2 = [0.2, 0.2, 0.1]; 
lightblue		= [0.08,	0.8,	1];
darkblue = [0.1, 0.01, 0.6]; darkblue2= [0.04, 0.06, 0.66];
verylight_blue     = [0.6, 0.9, 1];
lightgreen = [0.5, 0.98, 0.6]; lightgreen2= [0.45, 0.99, 0.55];
darkred = [0.9, 0.2, 0.1]; darkred2 = [0.9, 0.2, 0.1];

%%% figure 1:
 figure('Name', 'u and w functions');   

	plot(xaxis,u0, "-", 'Color',Blue, 'LineWidth', 2, 'MarkerSize',6); %% EE orig
	hold on;
  	% plot(xaxis,w0, "-", 'Color',Red, 'LineWidth', 2, 'MarkerSize',12); %% EE orig
	% hold on;	
    plot(xaxis,uex, ":p", 'Color',greenblue, 'LineWidth', 1.4, 'MarkerSize',7); %% EE orig
	hold on;	
    % plot(xaxis,wex, "--", 'Color',orange, 'LineWidth', 2.8, 'MarkerSize',12); %% EE orig
	% hold on;
    plot(xaxis,UC, '-s',  'Color',magenta, 'LineWidth', 2,'MarkerSize',6.5); %% Heun Ave
	hold on;
    plot(xaxis,ulh, '--',  'Color', black, 'LineWidth', 3.5, 'MarkerSize',9); %% Heun orig
 	hold on;
    
    XLabel = xlabel('x');
set(XLabel, 'FontSize', 18);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('u');
set(YLabel, 'FontSize', 18); 
set(YLabel, 'FontWeight', 'bold');

Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'u0','uexact','CCL','LH'},'Location','southeast','NumColumns', 2);
set(Legend, 'FontSize', 18); 
set(Legend, 'FontWeight', 'bold');
%title(Legend, 'My Legend Title');
% set(gca,'xScale','log');
% set(gca,'yScale','log');
grid on;

 	hold off;
    
%%% figure 2:
 figure('Name', 'Max Errors as a fuction of time step h');   

	plot(Axhstep(:),MaxD(:,1), ":s", 'Color',Red, 'LineWidth', 4.4, 'MarkerSize',11); %% EE orig
	hold on;
    plot(Axhstep(:),MaxD(:,3), '--h',  'Color', darkgreen, 'LineWidth', 3,'MarkerSize',9.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,5), '-.s',  'Color',darkblue, 'LineWidth', 3.4,'MarkerSize',5.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,7), '-o',  'Color', lightblue, 'LineWidth', 2.8, 'MarkerSize',15); %% Heun orig
 	hold on;
    plot(Axhstep(:),MaxD(:,9), ':*',  'Color', lightgreen, 'LineWidth', 4.3,'MarkerSize',9); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,11), '--<',  'Color', magenta, 'LineWidth', 2.6,'MarkerSize',7.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,13), '--m',  'Color', black, 'LineWidth', 3.4,'MarkerSize',4.5); %% Heun Ave
	hold on; 
    plot(Axhstep(:),MaxD(:,15), '-o',  'Color', darkred, 'LineWidth', 4,'MarkerSize',8.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,17), '-.h',  'Color',yellow, 'LineWidth', 4,'MarkerSize',6.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,19), '-x',  'Color', purple, 'LineWidth', 3, 'MarkerSize',15); %% Heun orig
 	hold on;
    plot(Axhstep(:),MaxD(:,2), '-d',  'Color', Red, 'LineWidth', 2.8, 'MarkerSize',13); %% EE AVE
	hold on;
    plot(Axhstep(:),MaxD(:,4), '--^',  'Color',  darkgreen2, 'LineWidth', 2.2,'MarkerSize',10); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,6), '-.X',  'Color', darkblue2, 'LineWidth', 2.9,'MarkerSize',10); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,8), '-',  'Color', verylight_blue, 'LineWidth', 4.2,'MarkerSize',5); %% Heun Ave
	hold on;
     plot(Axhstep(:),MaxD(:,10), '-^',  'Color', lightgreen2, 'LineWidth', 2.8,'MarkerSize',5.6); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,12), '-->',  'Color', darkredmagenta, 'LineWidth', 2.9,'MarkerSize',6.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,14), '-d',  'Color', black2, 'LineWidth', 2.9,'MarkerSize',8); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,16), '-O',  'Color', darkred2, 'LineWidth', 2.2,'MarkerSize',11); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,18), '--',  'Color', yellowdark, 'LineWidth', 4,'MarkerSize',7.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),MaxD(:,20), '-.p',  'Color', purple2, 'LineWidth', 2.4,'MarkerSize',7); %% Heun Ave
	hold on;
  	hold off;

XLabel = xlabel('Time step size \Delta\it{t}');
set(XLabel, 'FontSize', 18);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 18); 
set(YLabel, 'FontWeight', 'bold');

Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'FTCS: \it{u}','AB2: \it{u}','ADE: \it{u}','DF: \it{u}','OOEH: \it{u}','CpC: \it{u}','CCL: \it{u}','CLQ2: \it{u}','PI: \it{u}','LH: \it{u}','FTCS: \it{w}','AB2: \it{w}','ADE: \it{w}','DF: \it{w}','OOEH: \it{w}','CpC: \it{w}','CCL: \it{w}','CLQ2: \it{w}','PI: \it{w}','LH: \it{w}'},'Location','southeast','NumColumns', 2);
set(Legend, 'FontSize', 18); 
set(Legend, 'FontWeight', 'bold');
%title(Legend, 'My Legend Title');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;


%%%

%%% figure 3:
 figure('Name', 'L2D Errors as a fuction of time step h');   

	plot(Axhstep(:),L2D(:,1), ":s", 'Color',Red, 'LineWidth', 4.4, 'MarkerSize',11); %% EE orig
	hold on;
    plot(Axhstep(:),L2D(:,3), '--h',  'Color', darkgreen, 'LineWidth', 3,'MarkerSize',9.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,5), '-.s',  'Color',darkblue, 'LineWidth', 3.4,'MarkerSize',5.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,7), '-o',  'Color', lightblue, 'LineWidth', 2.8, 'MarkerSize',15); %% Heun orig
 	hold on;
    plot(Axhstep(:),L2D(:,9), ':*',  'Color', lightgreen, 'LineWidth', 4.3,'MarkerSize',9); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,11), '--<',  'Color', magenta, 'LineWidth', 2.6,'MarkerSize',7.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,13), '--m',  'Color', black, 'LineWidth', 3.4,'MarkerSize',4.5); %% Heun Ave
	hold on; 
    plot(Axhstep(:),L2D(:,15), '-o',  'Color', darkred, 'LineWidth', 4,'MarkerSize',8.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,17), '-.h',  'Color',yellow, 'LineWidth', 4,'MarkerSize',6.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,19), '-x',  'Color', purple, 'LineWidth', 3, 'MarkerSize',15); %% Heun orig
 	hold on;    
    plot(Axhstep(:),L2D(:,2), '-d',  'Color', Red, 'LineWidth', 2.8, 'MarkerSize',13); %% EE AVE
	hold on;
    plot(Axhstep(:),L2D(:,4), '--^',  'Color',  darkgreen2, 'LineWidth', 2.2,'MarkerSize',10); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,6), '-.X',  'Color', darkblue2, 'LineWidth', 2.9,'MarkerSize',10); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,8), '-',  'Color', verylight_blue, 'LineWidth', 4.2,'MarkerSize',5); %% Heun Ave
	hold on;
     plot(Axhstep(:),L2D(:,10), '-^',  'Color', lightgreen2, 'LineWidth', 2.8,'MarkerSize',5.6); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,12), '-->',  'Color', darkredmagenta, 'LineWidth', 2.9,'MarkerSize',6.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,14), '-d',  'Color', black2, 'LineWidth', 2.9,'MarkerSize',8); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,16), '-O',  'Color', darkred2, 'LineWidth', 2.2,'MarkerSize',11); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,18), '--',  'Color', yellowdark, 'LineWidth', 4,'MarkerSize',7.5); %% Heun Ave
	hold on;
    plot(Axhstep(:),L2D(:,20), '-.p',  'Color', purple2, 'LineWidth', 2.4,'MarkerSize',7); %% Heun Ave
	hold on;
  	hold off;

XLabel = xlabel('Time step size \Delta\it{t}');
set(XLabel, 'FontSize', 18);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 18); 
set(YLabel, 'FontWeight', 'bold');

Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'FTCS: \it{u}','AB2: \it{u}','ADE: \it{u}','DF: \it{u}','OOEH: \it{u}','CpC: \it{u}','CCL: \it{u}','CLQ2: \it{u}','PI: \it{u}','LH: \it{u}','FTCS: \it{w}','AB2: \it{w}','ADE: \it{w}','DF: \it{w}','OOEH: \it{w}','CpC: \it{w}','CCL: \it{w}','CLQ2: \it{w}','PI: \it{w}','LH: \it{w}'},'Location','southeast','NumColumns', 2);
set(Legend, 'FontSize', 18); 
set(Legend, 'FontWeight', 'bold');
%title(Legend, 'My Legend Title');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;




