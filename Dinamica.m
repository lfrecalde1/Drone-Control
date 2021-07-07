function DINAMICA_DRONE = Dinamica(Vref,v,PARAMETROS_DINAMICA,ts)
vref=[Vref(1);Vref(2);Vref(3);Vref(4)];
% matriz de inercia
x=PARAMETROS_DINAMICA;
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
    C12=x(8)+x(9)*v(4);
    C13=x(10);
    C14=x(11);
    C21=x(12)+x(13)*v(4);
    C22=x(14);
    C23=x(15);
    C24=x(16)+x(17)*v(4);
    C31=x(18);
    C32=x(19);
    C33=x(20);
    C34=x(21);
    C41=x(22);
    C42=x(23)+x(24)*v(4);
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
    vp = pinv(M)*(vref-C*v-G);
    v=vp*ts+v;
    DINAMICA_DRONE=v;
end
