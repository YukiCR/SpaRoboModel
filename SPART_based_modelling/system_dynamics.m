function ddx0 = system_dynamics(state,t, u)
L0=0.050;
L1=0.125;
L2=0.144;
L3=0.047;
L4=0.142;
L5=0.080;
L6=0.070;

%--- Manipulator Definition ----%
%Number of joints/links
data.n=5;

%First joint
data.man(1).type=1;
data.man(1).DH.d = L1;
data.man(1).DH.alpha = pi/2;
data.man(1).DH.a = 0;
data.man(1).DH.theta=0;
data.man(1).b = [0;L1/2;0];
data.man(1).mass=2;
data.man(1).I=diag([2,1,3])/10;

%Second joint
data.man(2).type=1;
data.man(2).DH.d = 0;
data.man(2).DH.alpha = 0;
data.man(2).DH.a = sqrt(L2^2+L3^2);
data.man(2).DH.theta=atan2(L2,L3);
data.man(2).b = [cos(-data.man(2).DH.theta),-sin(-data.man(2).DH.theta),0;sin(-data.man(2).DH.theta),cos(-data.man(2).DH.theta),0;0,0,1]*[L3^2/2;L2^2/2 + L3*L2;0]/(L2 + L3);
data.man(2).mass=2;
data.man(2).I=diag([2,1,3])/10;

%Third joint
data.man(3).type=1;
data.man(3).DH.d = 0;
data.man(3).DH.alpha = 0;
data.man(3).DH.a =L4;
data.man(3).DH.theta=-atan2(L2,L3);
data.man(3).b = [L4/2;0;0];
data.man(3).mass=2;
data.man(3).I=diag([2,1,3])/10;

%Fourth joint
data.man(4).type=1;
data.man(4).DH.d = 0;
data.man(4).DH.alpha = pi/2;
data.man(4).DH.a = 0;
data.man(4).DH.theta=pi/2;
data.man(4).b = [0;0;-L5/2];
data.man(4).mass=2;
data.man(4).I=diag([2,1,3])/10;

%Fifth joint
data.man(5).type=1;
data.man(5).DH.d = L5+L6;
data.man(5).DH.alpha =-pi/2;
data.man(5).DH.a = 0;
data.man(5).DH.theta=-pi/2;
data.man(5).b = [L6/2;0;0];
data.man(5).mass=2;
data.man(5).I=diag([2,1,3])/10;

%First joint location with respect to base
data.base.T_L0_J1=[eye(3),[0;0;L0];zeros(1,3),1];

%Base-spacecraft mass and inertia
data.base.mass=20;
data.base.I=diag([2,1,3]);

%End-Effector
data.EE.theta=-pi/2; %绕z轴再转-pi/2
data.EE.d=0;

q0 = state(1);
q1 = state(2);
q2 = state(3);
q3 = state(4);
x0 = state(5);
y0 = state(6);
z0 = state(7);
theta1 = state(8);
theta2 = state(9);
theta3 = state(10);
theta4 = state(11);
theta5 = state(12);
wx0 = state(13);
wy0 = state(14);
wz0 = state(15);
vx0 = state(16);
vy0 = state(17);
vz0 = state(18);
dtheta1 = state(19);
dtheta2 = state(20);
dtheta3 = state(21);
dtheta4 = state(22);
dtheta5 = state(23);

q = [q0;q1;q2;q3];
r0 = [x0;y0;z0];
R0 = rot_sys(q);  %Rotation from Base-spacecraft to inertial
w00 = [wx0;wy0;wz0];

%Joint variables [rad]
qm=[theta1;theta2;theta3;theta4;theta5];

%Velocities
u0=[wx0;wy0;wz0;vx0;vy0;vz0]; %Base-spacecraft velocity
um=[dtheta1;dtheta2;dtheta3;dtheta4;dtheta5]; %Joint velocities
dx0dt = [u0;um];

%--- Create robot structure ---%
[robot,T_Ln_EE] = DH_Serial2robot(data);

%--- Kinematics ---%
[RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,qm,robot);

%--- Differential Kinematics ---%
%Differential kinematics
[Bij,Bi0,P0,pm]=DiffKinematics(R0,r0,rL,e,g,robot);
%Velocities
[t0,tm]=Velocities(Bij,Bi0,P0,pm,u0,um,robot);

%--- Inertia Matrices ---%
%Inertias in inertial frames
[I0,Im]=I_I(R0,RL,robot);
%Mass Composite Body matrix
[M0_tilde,Mm_tilde]=MCB(I0,Im,Bij,Bi0,robot);
%Generalized Inertia matrix
[H0, H0m, Hm] = GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
%Generalized Convective Inertia matrix
[C0, C0m, Cm0, Cm] = CIM(t0,tm,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);

H = [H0,H0m;H0m',Hm];
C = [C0,C0m;Cm0,Cm];
d = [0;0;0;0;0;0;0;0;0;0;0];
ddx0dt = inv(H)*(d+u-C*dx0dt);

ddq1 = ddx0dt(1);
ddq2 = ddx0dt(2);
ddq3 = ddx0dt(3);
ddq4 = ddx0dt(4);
ddq5 = ddx0dt(5);
ddq6 = ddx0dt(6);
ddq7 = ddx0dt(7);
ddq8 = ddx0dt(8);
ddq9 = ddx0dt(9);
ddq10 = ddx0dt(10);
ddq11 = ddx0dt(11);

qv = [q1; q2; q3];
dq = 0.5 * [-qv'; (q0 * eye(3) - ccsz(qv))'] * w00;
dq0dt = dq(1);
dq1dt = dq(2);
dq2dt = dq(3);
dq3dt = dq(4);
ddx0 = [dq0dt;dq1dt;dq2dt;dq3dt;vx0;vy0;vz0;dtheta1;dtheta2;dtheta3;dtheta4;dtheta5;ddq1;ddq2;ddq3;ddq4;ddq5;ddq6;ddq7;ddq8;ddq9;ddq10;ddq11];

end
