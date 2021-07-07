%% CONTROLADOR BASADO EN OPTIMIZACION DE UN DRONE

%% LIMPIAR VARIABLES DEL SISTEMA
clc,clear all,close all;

%% DEFINICION DEL TIEMPO DE SIMULACION
to=0;
t_final=100;
ts=0.1;
t=[to:ts:t_final];
%% DEFINICION DE LAS POSICIONES INICIALES DEL ROBOT AEREO

hx(1)=1;
hy(1)=1;
hz(1)=1;
hth(1)=45*pi/180;
%% DEFINICION DE LAS VELOCIDADES INICIALES DEL ROBOT AEREO


hxp(1)=0;
hyp(1)=0;
hzp(1)=0;
hthp(1)=0;

%% DEFINICION DE LAS SENALES DE CONTROL DEL SISTEMA


%% DEFICNION DE LAS VELOCIDADES DE CONTROL REALES 
ul(1)=0;
um(1)=0;
un(1)=0;
w(1)=0;

%% SENALES DESEADAS DEL SISTEMA
hxd=0.25*t+2;
hxdp=0.25*ones(1,length(t));
hxdpp=0*ones(1,length(t));

hyd=2*sin(t/8)+0.05*t-4;
hydp=(1/8)*2*cos(t/8)+0.05;
hydpp=-(1/8)*(1/8)*2*sin(t/8);

hzd=10+1.5*sin(t/10);
hzdp=1.5*(1/10)*cos(t/10);
hzdpp=-1.5*(1/10)*1/10*sin(t/10);

hthd=(atan2(hydp,hxdp));
hthdp=diff([0 hthd])/ts;

%% GENERACION DE LAS ACCIONES DE CONTROL DE REFERENCIA
vRef = sqrt(hxdp.^2+hydp.^2);
wRef = (hxdp.*hydpp-hydp.*hxdpp)./(hxdp.^2+hydp.^2);

%% DEFICNION DE LAS MATRICES PARA EL CONTROL CINEMATICO
K1=diag([1 1 1 1]);
K2=diag([1 1 1 1]);

%% CARGAR VALORES DE LA DINAMICA DEL DRONE
load('DINAMICA_DRONE.mat')
chi_real=x'*ones(1,length(t));
chi_estimado(:,1)=x;

for k=1:length(t)-4
    tic;
    
    %% VECTORES DE LOS ERRORES PRIMARIOS
    hxe(k)=hxd(k)-hx(k);
    hye(k)=hyd(k)-hy(k);
    hze(k)=hzd(k)-hz(k);
    hthe(k)=wrapToPi(hthd(k)-hth(k));
 
    %% VECTORES DE ESTADOS EN FORMA GENERAL CONTROL PRIMARIO
    h=[hx(k) hy(k) hz(k) hth(k)]';
    hd=[hxd(k) hyd(k) hzd(k) hthd(k)]';
    he=hd-h;
    hdp=[hxdp(k) hydp(k) hzdp(k) hthdp(k)]';

   %% DEFINICION DEL JACOBIANO DEL SISTEMA
   J11=cos(hth(k));
   J12=-sin(hth(k));
   J13=0;
   J14=0;
   
   J21=sin(hth(k));
   J22=cos(hth(k));
   J23=0;
   J24=0;
   
   J31=0;
   J32=0;
   J33=1;
   J34=0;
   
   J41=0;
   J42=0;
   J43=0;
   J44=1;
   
   J=[J11,J12,J13,J14;...
      J21,J22,J23,J24;...
      J31,J32,J33,J34;...
      J41,J42,J43,J44]; 
  qref_c=pinv(J)*(hdp+K2*tanh(pinv(K2)*K1*he));
  
  ulref_c(k)=qref_c(1);
  umref_c(k)=qref_c(2);
  unref_c(k)=qref_c(3);
  wref_c(k)=qref_c(4);
  
  %% GENERACION DE LOS VECTORES PARA LA DINAMICA Y COMPENSACION DINAMICA
  %% ERRORES DE VELOCIDAD VREF
  ule(k)=ulref_c(k)-ul(k);
  ume(k)=umref_c(k)-um(k);
  une(k)=unref_c(k)-un(k);
  we(k)=wref_c(k)-w(k);
  
  vref_e=[ule(k) ume(k) une(k) we(k)]';
  
  %% DERIVADAS DE LAS VREF DESEADAS
  vrefp_ul=diff([ulref_c ulref_c(end)])/ts;
  vrefp_um=diff([umref_c umref_c(end)])/ts;
  vrefp_un=diff([unref_c unref_c(end)])/ts;
  vrefp_w=diff([wref_c wref_c(end)])/ts;

  vrefp=[vrefp_ul(k) vrefp_um(k) vrefp_un(k) vrefp_w(k)]';
  
  v=[ul(k);um(k);un(k);w(k)];
  
  %% COMPENSACION DINAMICA PARA EL DRONE
 COMPENSACION=Compensacion(vrefp,vref_e,v,chi_real(:,k),chi_estimado(:,k),ts);
 
 ulref(k)=COMPENSACION(1);
 umref(k)=COMPENSACION(2);
 unref(k)=COMPENSACION(3);
 wref(k)=COMPENSACION(4);
 chi_estimado(:,k+1)=COMPENSACION(5:end,1);

 vref=[ulref(k);umref(k);unref(k);wref(k)];
 %% DINAMICA DEL DORNE REAL ESTIMADA
 DINAMICA_DRONE = Dinamica(vref,v,chi_real(:,k),ts);
  
 ul(k+1)=DINAMICA_DRONE(1)+rand(1)*0.1;
 um(k+1)=DINAMICA_DRONE(2)+rand(1)*0.1;
 un(k+1)=DINAMICA_DRONE(3)+rand(1)*0.1;
 w(k+1)= DINAMICA_DRONE(4)+rand(1)*0.1;
 
 hth(k+1)=wrapToPi(hth(k)+w(k+1)*ts);
 
  qreal=[ul(k+1);um(k+1);un(k+1);w(k+1)];
  %% EVOLUCION DEL SISTEMA
  hxp(k+1)=cos(hth(k+1))*ul(k+1)-sin(hth(k+1))*um(k+1);
  hyp(k+1)=sin(hth(k+1))*ul(k+1)+cos(hth(k+1))*um(k+1);
  hzp(k+1)=un(k+1);
 
  %% INTEGRACION NUMERICA PARA CONOCER LOS ESTADOS DEL SISTEMA
  hx(k+1)=hx(k)+hxp(k+1)*ts;
  hy(k+1)=hy(k)+hyp(k+1)*ts;
  hz(k+1)=hz(k)+hzp(k+1)*ts;

 %% ALMACENAMIENTO DEL TIEMPO DE SAMPLEO
  t_sample(k)=toc;
end
close all; paso=1; 
%a) Parámetros del cuadro de animación
figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 3]);
h = light;
h.Color=[0.65,0.65,0.65];
h.Style = 'infinite';
%b) Dimenciones del Robot
   Drone_Parameters(0.02);
%c) Dibujo del Robot    
    G2=Drone_Plot_3D(hx(1),hy(1),hz(1),0,0,hth(1));hold on

    plot3(hx(1),hy(1),hz(1),'--','Color',[56,171,217]/255,'linewidth',1.5);hold on,grid on   
    plot3(hxd(1),hyd(1),hzd(1),'Color',[32,185,29]/255,'linewidth',1.5);

view(20,15);
for k = 1:20:length(t)-4
    drawnow
    delete(G2);
   
    G2=Drone_Plot_3D(hx(k),hy(k),hz(k),0,0,hth(k));hold on
    
    plot3(hxd(1:k),hyd(1:k),hzd(1:k),'Color',[32,185,29]/255,'linewidth',1.5);
    plot3(hx(1:k),hy(1:k),hz(1:k),'--','Color',[56,171,217]/255,'linewidth',1.5);
    
    legend('$\mathbf{h}$','$\mathbf{h}_{des}$','Interpreter','latex','FontSize',11,'Location','northwest','Orientation','horizontal');
    legend('boxoff')
    title('$\textrm{Movement Executed by the Aerial Robot}$','Interpreter','latex','FontSize',11);
    xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);
    
end
print -dpng SIMULATION_1
print -depsc SIMULATION_1

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
subplot(2,1,1)
plot(t(1:length(hxe)),hxe,'Color',[226,76,44]/255,'linewidth',1); hold on;
plot(t(1:length(hxe)),hye,'Color',[46,188,89]/255,'linewidth',1); hold on;
plot(t(1:length(hxe)),hze,'Color',[26,115,160]/255,'linewidth',1);hold on;
plot(t(1:length(hxe)),hthe,'Color',[83,57,217]/255,'linewidth',1);hold on;
grid on;
legend('$\tilde{h_{x}}$','$\tilde{h_{y}}$','$\tilde{h_{z}}$','$\tilde{h_{\psi}}$','Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Evolution of Control Errors}$','Interpreter','latex','FontSize',9);
ylabel('$[m]$','Interpreter','latex','FontSize',9);
subplot(2,1,2)
plot(t(1:length(ule)),ule,'Color',[226,76,44]/255,'linewidth',1); hold on;
plot(t(1:length(ule)),ume,'Color',[46,188,89]/255,'linewidth',1); hold on;
plot(t(1:length(ule)),une,'Color',[26,115,160]/255,'linewidth',1);hold on;
plot(t(1:length(ule)),we,'Color',[83,57,217]/255,'linewidth',1);hold on;
grid on;
legend('$\tilde{\mu_{l}}$','$\tilde{\mu_{m}}$','$\tilde{\mu_{n}}$','$\tilde{\omega}$','Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Evolution of Control Errors}$','Interpreter','latex','FontSize',9);
ylabel('$[m]$','Interpreter','latex','FontSize',9);
print -dpng ERROR_SISTEMA
print -depsc ERROR_SISTEMA

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
subplot(2,1,1)
plot(t(1:length(ulref_c)),ulref_c,'Color',[226,76,44]/255,'linewidth',1); hold on
plot(t(1:length(ulref_c)),umref_c,'Color',[46,188,89]/255,'linewidth',1); hold on
plot(t(1:length(ulref_c)),unref_c,'Color',[26,115,160]/255,'linewidth',1); hold on
plot(t(1:length(ulref_c)),wref_c,'Color',[83,57,217]/255,'linewidth',1); hold on
plot(t(1:length(ul)),ul,'--','Color',[226,76,44]/255,'linewidth',1); hold on
plot(t(1:length(ul)),um,'--','Color',[46,188,89]/255,'linewidth',1); hold on
plot(t(1:length(ul)),un,'--','Color',[26,115,160]/255,'linewidth',1); hold on
plot(t(1:length(ul)),w,'--','Color',[83,57,217]/255,'linewidth',1); hold on
grid on;
legend('$\mu_{lc}$','$\mu_{mc}$','$\mu_{nc}$','$\omega_{c}$','$\mu_{l}$','$\mu_{m}$','$\mu_{n}$','$\omega$','Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Control Values}$','Interpreter','latex','FontSize',9);
ylabel('$[rad/s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);
subplot(2,1,2)
plot(t(1:length(ulref)),ulref,'Color',[226,76,44]/255,'linewidth',1); hold on
plot(t(1:length(ulref)),umref,'Color',[46,188,89]/255,'linewidth',1); hold on
plot(t(1:length(ulref)),unref,'Color',[26,115,160]/255,'linewidth',1); hold on
plot(t(1:length(ulref)),wref,'Color',[83,57,217]/255,'linewidth',1); hold on
grid on;
legend('$\mu_{lref}$','$\mu_{mref}$','$\mu_{nref}$','$\omega_{ref}$','Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Control Values}$','Interpreter','latex','FontSize',9);
ylabel('$[rad/s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);
print -dpng CONTROL_VALUES
print -depsc CONTROL_VALUES

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(t_sample)),t_sample,'Color',[46,188,89]/255,'linewidth',1); hold on
grid on;
legend('$t_{sample}$','Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Control Values}$','Interpreter','latex','FontSize',9);
ylabel('$[s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);  
print -dpng SAMPLE_TIME
print -depsc SAMPLE_TIME
figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(chi_estimado)),chi_estimado,'linewidth',1);
hold on
grid on
title('$\textrm{Adaptative Parameters}$','Interpreter','latex','FontSize',9);
legend('$\chi_e$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Parameters}$','Interpreter','latex','FontSize',9);
print -dpng PARAMETROS
print -depsc PARAMETROS