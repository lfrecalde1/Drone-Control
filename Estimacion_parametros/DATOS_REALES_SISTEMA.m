clc,clear all,close all;
warning off
global t  UL UM UN W ULP UMP UNP WP WREF ULREF UMREF UNREF
load('DATOS_REALES.mat')

ul_sal=xu_p;
um_sal=yu_p;
un_sal=zu_p;
w_sal=w;

ULREF=ul_ref(1:length(ul_sal));
UMREF=um_ref(1:length(ul_sal));
UNREF=un_ref(1:length(ul_sal));
WREF=w_ref(1:length(ul_sal));

t=0:0.1:(length(ULREF)-1)*0.1;

UL=ul_sal;
UM=um_sal;
UN=un_sal;
W=w_sal;

ULP=[0 diff(UL)/ts];
UMP=[0 diff(UM)/ts];
UNP=[0 diff(UN)/ts];
WP=[0 diff(W)/ts];

landa=1;%lambda
F1=tf(landa,[1 landa]);


UL=lsim(F1,UL,t);
UM=lsim(F1,UM,t);
UN=lsim(F1,UN,t);
W=lsim(F1,W,t);

ULP=lsim(F1,ULP,t);
UMP=lsim(F1,UMP,t);
UNP=lsim(F1,UNP,t);
WP=lsim(F1,WP,t);

ULREF=lsim(F1,ULREF,t);
UMREF=lsim(F1,UMREF,t);
UNREF=lsim(F1,UNREF,t);
WREF=lsim(F1,WREF,t);

x0=zeros(1,27)+0.1;

options = optimoptions('fminunc','Display','iter',...
    'Algorithm','quasi-newton','TolFun',1e-20,'TolX',1e-20,...
    'MaxFunEvals',1e+100,'MaxIter',1000);
[x,fval,exit,salida] = fminunc(@Dinamica_estimacion,x0,options);

x;
v(1,1)=0;
v(2,1)=0;
v(3,1)=0;
v(4,1)=0;

va(1,1)=0;
va(2,1)=0;
va(3,1)=0;
va(4,1)=0;

save('DINAMICA_DRONE.mat','x')

for i=1:length(ULREF)
    v=[ul_sal(i);um_sal(i);un_sal(i);w_sal(i)];
    vref=[ul_ref(i);um_ref(i);un_ref(i);w_ref(i)];
    
 
    % matriz de inercia
    M11=x(1);
    M12=0;
    M13=0;
    M14=x(2);
    M21=0;
    M22=x(3);
    M23=0;
    M24=0;
    M31=0;
    M32=0;
    M33=x(4);
    M34=0;
    M41=x(5);
    M42=0;
    M43=0;
    M44=x(6);
    
    
   
    M=[M11,M12,M13,M14;...
       M21,M22,M23,M24;...
       M31,M32,M33,M34;...
       M41,M42,M43,M44];
    % matriz centifuga-centripeta
    C11=x(7);
    C12=x(8)+x(9)*va(4,i);
    C13=x(10);
    C14=x(11);
    C21=x(12)+x(13)*va(4,i);
    C22=x(14);
    C23=x(15);
    C24=x(16)+x(17)*va(4,i);
    C31=x(18);
    C32=x(19);
    C33=x(20);
    C34=x(21);
    C41=x(22);
    C42=x(23)+x(24)*va(4,i);
    C43=x(25);
    C44=x(26);

   
    
     C=[C11,C12,C13,C14;...
       C21,C22,C23,C24;...
       C31,C32,C33,C34;...
       C41,C42,C43,C44];
    % vector de gravedad
    
G11=0;
G21=0;
G31=x(27);
G41=0;
G=[G11;G21;G31;G41];
    %Dinamica
    vp = pinv(M)*(vref-C*va(:,i)-G);
    va(:,i+1)=vp*ts+va(:,i);
    
    end
figure(1)
% Sentido Ley de la mano derech
subplot(4,1,1)
plot(ul_ref,'g') % ul(k)
hold on
grid on
plot(va(1,1:length(ULREF)),'b')  % ul_real(k)
plot(ul_sal,'r')
legend("ul ref","ul estimado",'ul real')
ylabel('x [m/s]'); xlabel('time [ms]');
title ("VELOCIDADES DE MANIPULABILIDAD ul, um, un, w")

subplot(4,1,2)
plot(um_ref,'g')    % um(k)
hold on
grid on
plot(va(2,1:length(UMREF)),'b') % um_real(k)  
plot(um_sal,'r')
legend("um ref","um estimado",'um real')
ylabel('y [m/s]'); xlabel('time [ms]')

subplot(4,1,3)
plot(un_ref,'g')    % um(k)
hold on
grid on
plot(va(3,1:length(UNREF)),'b') % um_real(k)  
plot(un_sal,'r')
legend("un ref","un estimado",'un real')
ylabel('z [m/s]'); xlabel('time [ms]')

subplot(4,1,4)
plot(w_ref,'g')    % um(k)
hold on
grid on
plot(va(4,1:length(UNREF)),'b') % um_real(k)  
plot(w_sal,'r')
legend("w ref","w estimado",'w real')
ylabel('omega [rad/s]'); xlabel('time [ms]')
print -dpng MODELACION_DINAMICA
print -depsc MODELACION_DINAMICA