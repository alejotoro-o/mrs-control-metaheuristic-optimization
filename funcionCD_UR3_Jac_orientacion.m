%% PRUEBA
% clear
% clc
% 
% q=deg2rad([0, -90, -90, -90, 90,0]); %para prueba 1
% % q=deg2rad([0.000061, -95.628528, -56.970721, -117.400972, 89.999779, 0.000061]); %Prueba 2
% % q=[0.281207920409350	-1.793650390169974	-1.393057747583438	0.699452335409308	0.608085339431503	2.003472748725389];
% q1=q(1); q2=q(2); q3=q(3); q4=q(4);q5=q(5);q6=q(6);
% error=funcionCD_UR3(q1,q2,q3,q4,q5,q6)
%%
function [error_f]=funcionCD_UR3_Jac_orientacion(var)

%Parámetros a estimar
q1=var(1);
q2=var(2);
q3=var(3);
q4=var(4);
q5=var(5);
q6=var(6);

% -----------------------------------------------------------
%           POSE DESEADA - Caso 1: "HOME" de pick and place
% -----------------------------------------------------------
xd=0.295400;        yd= -0.110400;     zd= 0.313150; 
rxd=deg2rad(127.279);        ryd= deg2rad(-127.279);     rzd=0; %orientación para el robot EF "viendo hacia abajo"

% -----------------------------------------------------------
%          POSE DESEADA - Caso 2 "Home" de pintura
% -----------------------------------------------------------
% xd=0.295400;        yd= -0.110400;     zd= 0.478950; 
% rxd=deg2rad(-69.282);        ryd= deg2rad(69.282);     rzd=-69.282;

% -----------------------------------------------------------
%          POSE DESEADA - Caso 3: Cercano a singularidad
% -----------------------------------------------------------
% xd=0.295400;        yd= -0.110400;     zd= 0.638;
% rxd=deg2rad(-69.282);        ryd= deg2rad(69.282);     rzd=-69.282; %orientación para el robot EF "viendo hacia abajo"


% -----------------------------------------------------------
%                  PARÁMETROS DH - UR3
% -----------------------------------------------------------

d1=0.1519;      d2=0;              d3=0;
d4=0.11235;     d5=0.08535;         d6=0.0819;
     
alp1=pi/2;      alp2=0;             alp3=0;
alp4=pi/2;      alp5=-pi/2;         alp6=0;
    
a1=0;           a2=-0.24365;        a3=-0.21325;
a4=0;           a5=0;               a6=0;
% -----------------------------------------------------------
%                  CD - UR3
% -----------------------------------------------------------

dh1=[cos(q1)	-sin(q1)*cos(alp1)	sin(q1)*sin(alp1)	a1*cos(q1)
    sin(q1)  	cos(q1)*cos(alp1)  	-cos(q1)*sin(alp1)	a1*sin(q1)
    0	sin(alp1)	cos(alp1)	d1
    0	0	0	1];

dh2=[cos(q2)	-sin(q2)*cos(alp2)	sin(q2)*sin(alp2)	a2*cos(q2)
    sin(q2)  	cos(q2)*cos(alp2)  	-cos(q2)*sin(alp2)	a2*sin(q2)
    0	sin(alp2)	cos(alp2)	d2
    0	0	0	1];

dh3=[cos(q3)	-sin(q3)*cos(alp3)	sin(q3)*sin(alp3)	a3*cos(q3)
    sin(q3)  	cos(q3)*cos(alp3)  	-cos(q3)*sin(alp3)	a3*sin(q3)
    0	sin(alp3)	cos(alp3)	d3
    0	0	0	1];

dh4=[cos(q4)	-sin(q4)*cos(alp4)	sin(q4)*sin(alp4)	a4*cos(q4)
    sin(q4)  	cos(q4)*cos(alp4)  	-cos(q4)*sin(alp4)	a4*sin(q4)
    0	sin(alp4)	cos(alp4)	d4
    0	0	0	1];

dh5=[cos(q5)	-sin(q5)*cos(alp5)	sin(q5)*sin(alp5)	a5*cos(q5)
    sin(q5)  	cos(q5)*cos(alp5)  	-cos(q5)*sin(alp5)	a5*sin(q5)
    0	sin(alp5)	cos(alp5)	d5
    0	0	0	1];

dh6=[cos(q6)	-sin(q6)*cos(alp6)	sin(q6)*sin(alp6)	a6*cos(q6)
    sin(q6)  	cos(q6)*cos(alp6)  	-cos(q6)*sin(alp6)	a6*sin(q6)
    0	sin(alp6)	cos(alp6)	d6
    0	0	0	1];

DHF=dh1*dh2*dh3*dh4*dh5*dh6;
% -----------------------------------------------------------
%                  POSICIÓN EF - UR3
% -----------------------------------------------------------
xrobot=DHF(1,4);
yrobot=DHF(2,4);
zrobot=DHF(3,4);
% -----------------------------------------------------------
%                  ORIENTACIÓN EF - UR3
% -----------------------------------------------------------
R = DHF(1:3, 1:3);    
ur_vector_rad = rotationToURVector(R);
ur_vector_deg = rad2deg(ur_vector_rad);
orientacion=ur_vector_deg;
rx=orientacion(1);
ry=orientacion(2);
rz=orientacion(3);

% -----------------------------------------------------------
%                  POSE ACTUAL - UR3
% -----------------------------------------------------------

X=[xrobot,yrobot,zrobot,orientacion];

% -----------------------------------------------------------
%                  ÇÁLCULO DE ERRORES CARTESIANOS- UR3
% -----------------------------------------------------------

% errores de posicion
ex=xd-xrobot;   ey=yd-yrobot;     ez=zd-zrobot;
ep=[ex;ey;ez];

% -----------------------------------------------------------
%                   CÁLCULO DE JACOBIANO
% -----------------------------------------------------------

% Transformaciones intermedias desde la base
T00 = eye(4);
T01 = dh1;
T02 = T01 * dh2;
T03 = T02 * dh3;
T04 = T03 * dh4;
T05 = T04 * dh5;
T06 = T05 * dh6; % Es lo mismo que DHF

% Extracción de vectores de posicion (o) y orientacion (z) de cada S.Ref
o0 = T00(1:3, 4);
o1 = T01(1:3, 4);
o2 = T02(1:3, 4);
o3 = T03(1:3, 4);
o4 = T04(1:3, 4);
o5 = T05(1:3, 4);
o_end = T06(1:3, 4);

z0 = T00(1:3, 3);
z1 = T01(1:3, 3);
z2 = T02(1:3, 3);
z3 = T03(1:3, 3);
z4 = T04(1:3, 3);
z5 = T05(1:3, 3);

% Calcular cada columna del jacobiano
J1 = [cross(z0, o_end - o0); z0];
J2 = [cross(z1, o_end - o1); z1];
J3 = [cross(z2, o_end - o2); z2];
J4 = [cross(z3, o_end - o3); z3];
J5 = [cross(z4, o_end - o4); z4];
J6 = [cross(z5, o_end - o5); z5];

% Ensamblar el jacobiano final
J = [J1, J2, J3, J4, J5, J6];

% -----------------------------------------------------------
%           ÇÁLCULO DE ERRORES ORIENTACIÓN - UR3
% -----------------------------------------------------------

rd=[rxd,ryd,rzd];
Rd=URVectorToRotation(rd);
% errores de orientacion
% R: actual, R_d: deseada
R_e = Rd.' * R;
ew = 0.5 * [R_e(3,2)-R_e(2,3);
             R_e(1,3)-R_e(3,1);
             R_e(2,1)-R_e(1,2)];

% -----------------------------------------------------------
%           ÇÁLCULO DE ERROR DE POSE - UR3
% -----------------------------------------------------------
% ep está en metros.
% ew está en radianes.
% Al sumarlos así se mezclan unidades. 
% para que tenga significado físico, puedes ponderar la parte angular con un factor de escala 
% λ (m/rad). λ=0.1 m/rad → 1 rad de error angular ≈ 10 cm de error lineal.
%            λ=0.05 m/rad → 1 rad de error angular ≈ 5 cm.
% Esto hace que un error razonable en orientación (ej: 0.1 rad = ~5.7°) se “vea” como 0.01 m de error equivalente.

lambda = 0.1;   % metros por radian
e_s = [ep; lambda*ew];

% -----------------------------------------------------------
%      ÇÁLCULO DE NORMA DE ERROR CON SINGULARIDADES - UR3
% -----------------------------------------------------------

Jdet=det(J); %Calculo det del Jacobiano, si det(J)=0 es singularidad
if Jdet==0
    error_f=norm(e_s)+1000;
else
    error_f=norm(e_s);
end

% -----------------------------------------------------------
%      FUNCIONES USADAS
% -----------------------------------------------------------

%Para crear vector de rotación deseado a Matriz Rd 3x3
function R = URVectorToRotation(r)
    theta = norm(r);
    if theta < 1e-12
        R = eye(3);
        return;
    end
    k = r / theta;  % eje unitario
    K = [  0   -k(3)  k(2);
          k(3)   0   -k(1);
         -k(2)  k(1)   0 ];
    R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);
end

%Para ángulos en formato UR
function ur_rotation_vector = rotationToURVector(R)
    % Verificar que la matriz sea válida
    if ~isequal(size(R), [3, 3])
        error('La matriz de rotación debe ser 3x3.');
    end
    
    % Convertir matriz de rotación a ángulo-eje
    axang = rotm2axang(R);
    axis = axang(1:3);  % Eje de rotación
    angle = axang(4);   % Ángulo de rotación
    
    % Escalar el eje por el ángulo
    ur_rotation_vector = axis * angle;
end

end