machinery = readtable('+user/3bus_PARK/machinery.xlsx');

parameter = machinery(:, {'Xd', 'Xd_prime', 'Xq', 'Tdo', 'M', 'D'});
parameter = parameter.Variables;

Xd = parameter(1,1);
Xdp = parameter(1,2);
Xq = parameter(1,3);
Tdo = parameter(1,4);
M = parameter(1,5);
D = parameter(1,6);

%AVR
Ka = 20; Te = 0.05;

%internal
omega0 = 2*pi*60;
delta = net.x_equilibrium(1);
E = net.x_equilibrium(3);
Vabs = abs(net.V_equilibrium(1));
Vangle = angle(net.V_equilibrium(1));
Iabs = abs(net.I_equilibrium(1));
Iangle = angle(net.I_equilibrium(1));

%Define syms
syms s_delta s_omega s_E s_Vfd s_Vpss s_Pmech
%Define Vd, Vq, Id, Iq
Vq = Vabs*cos(s_delta-Vangle);
Vd = Vabs*sin(s_delta-Vangle);
%Iq = Iabs*cos(s_delta-Iangle);
%Id = Iabs*sin(s_delta-Iangle);

Id = (s_E-Vq)/Xdp;
Iq = Vd/Xq;

%PSS
[A_pss,B_pss,C_pss,D_pss] = ssdata(net.a_controller_local{1}.pss.get_sys);

%1axis, num_x = 7

A = [0 omega0 0 0 0 0 0 ;
     0 -D/M   0 0 0 0 0;
     0  0     0 0 0 0 0;
     0  D_pss*(Ka/Te)  0 -1/Te C_pss.*(Ka/Te);
     [zeros(3,1),B_pss,zeros(3,2),A_pss]];

L = [0 0 0 0 0 0;
     -1/M 0 0 0 0 1/M;
     0 -1/Tdo 0 1/Tdo 0 0;
     0 0 -Ka/Te 0 -Ka/Te 0;zeros(3,6);
     ];

gamma = [eye(4,4),zeros(4,3)];


%too lazy to differentiate P with delta and E
%{
temp = (1/Xdp - 1/Xq)*Vabs^2*((sin(delta-Vangle))^2-(cos(delta-Vangle))^2);
environ = [E*Vabs*cos(delta-Vangle)/Xd+temp 0 Vabs*sin(delta-Vangle)/Xdp 0;
           (Xd/Xdp-1)*Vabs*sin(delta-Vangle) 0 Xd/Xdp 0;
           zeros(1,4);
           0 0 0 1;
           zeros(2,4)];
%}

Efd = Xd*s_E/Xdp - (Xd/Xdp-1)*Vq;
P = Vq*Iq+Vd*Id;
Environ = [ P; Efd; Vabs; s_Vfd; s_Vpss; s_Pmech];
J = jacobian(Environ, [s_delta, s_omega, s_E, s_Vfd]);
environ = eval(subs(J,[s_delta,s_E],[delta,E]));

A_ = (L*environ*gamma);
Ahat = A+A_;
Bw = L*environ;

%2axis, num_x = 8
%{
A = [0 omega0 0 0 0 0 0 0;
     0 -D/M   0 0 0 0 0 0;
     0  0     0 0 0 0 0 0;
     0  0  0 -1/Te 0 0 0;
     0  0  0   0   0 0 0;
     0  0  0   0   0 0 0;
     0  0  0   0   0 0 0;];
%}



