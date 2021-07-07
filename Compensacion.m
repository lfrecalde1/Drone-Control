function COMPENSACION= Compensacion(vrefp,vref_e,v,PARAMETROS_DINAMICA,PARAMETROS_ESTIMADOS,ts)
% matriz de inercia
%% DATOS DE LAS VELOCIDADES DEL DRONE
mu_l=v(1);
mu_m=v(2);
mu_n=v(3);
omega=v(4);
F=PARAMETROS_DINAMICA;
%% PARAMETROS PARA LA ADAPTACION DE LA DINAMICA
chi_real=PARAMETROS_DINAMICA;
chi_estimado=PARAMETROS_ESTIMADOS;
%% ERROR DE PARAMETROS D ELA DINAMICA
 error_chi=chi_estimado-chi_real;
 %% GANACIAS DE LA ADAPTACION DE PARAMETROS
  CHI=1*eye(27);
  H=1*eye(4);
  T=1*eye(27);
  
    M11=F(1);
    M12=0;
    M13=0;
    M14=F(2);
    M21=0;
    M22=F(3);
    M23=0;
    M24=0;
    M31=0;
    M32=0;
    M33=F(4);
    M34=0;
    M41=F(5);
    M42=0;
    M43=0;
    M44=F(6);
    
    
   
    M=[M11,M12,M13,M14;...
       M21,M22,M23,M24;...
       M31,M32,M33,M34;...
       M41,M42,M43,M44];
    % matriz centifuga-centripeta
    C11=F(7);
    C12=F(8)+F(9)*v(4);
    C13=F(10);
    C14=F(11);
    C21=F(12)+F(13)*v(4);
    C22=F(14);
    C23=F(15);
    C24=F(16)+F(17)*v(4);
    C31=F(18);
    C32=F(19);
    C33=F(20);
    C34=F(21);
    C41=F(22);
    C42=F(23)+F(24)*v(4);
    C43=F(25);
    C44=F(26);

   
    
     C=[C11,C12,C13,C14;...
       C21,C22,C23,C24;...
       C31,C32,C33,C34;...
       C41,C42,C43,C44];
    % vector de gravedad
G11=0;
G21=0;
G31=F(27);
G41=0;
G=[G11;G21;G31;G41];
%% GANANCIAS DL CONTROLADOR DINAMICO
K1=diag([2 2 2 2]);
K2=diag([1 1 1 1]);

%% COMPENSACION DINAMICA
control=vrefp+K2*tanh(inv(K2)*K1*vref_e);
mu_lp=control(1);
mu_mp=control(2);
mu_np=control(3);
omega_p=control(4);
    
n=[mu_lp, omega_p, 0, 0, 0, 0, mu_l, mu_m, mu_m*omega, mu_n, omega, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
   0, 0, mu_mp, 0, 0, 0, 0, 0, 0, 0, 0, mu_l, mu_l*omega, mu_m, mu_n, omega, omega^2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
   0, 0, 0, mu_np, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, mu_l, mu_m, mu_n, omega, 0, 0, 0, 0, 0, 1;...
   0, 0, 0, 0, mu_lp, omega_p, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, mu_l, mu_m, mu_m*omega, mu_n, omega, 0];
   
DI= n*chi_real+n*error_chi;

chi_p=inv(CHI)*n'*H*vref_e-inv(CHI)*T*error_chi;
chi_estimado=chi_estimado+chi_p*ts;

COMPENSACION=[DI;chi_estimado];
end

