machinery = readtable('+user/3bus_PARK/machinery.xlsx');

parameter = machinery(:, {'Xd', 'Xd_prime', 'Xd_pp','Xq', 'Xq_prime','Xq_pp','Tdo', 'Tqo','TTdo','TTqo','X_ls','M', 'D'});
parameter = parameter.Variables;

Xd = parameter(1,1);
Xdp = parameter(1,2);
Xdpp = parameter(1,3);
Xq = parameter(1,4);
Xqp = parameter(1,5);
Xqpp = parameter(1,6);
Tdo = parameter(1,7);
Tqo = parameter(1,8);
TTdo = parameter(1,9);
TTqo = parameter(1,10);
Xls = parameter(1,11);
M = parameter(1,12);
D = parameter(1,13);

%AVR
Ka = 20; Te = 0.05;

%internal
omega0 = 2*pi*60;
delta = net.x_equilibrium(1);
Eq = net.x_equilibrium(3);
Ed = net.x_equilibrium(4);
psiq = net.x_equilibrium(5);
psid = net.x_equilibrium(6);
Vabs = abs(net.V_equilibrium(1));
Vangle = angle(net.V_equilibrium(1));
Iabs = abs(net.I_equilibrium(1));
Iangle = angle(net.I_equilibrium(1));

%PSS
[A_pss,B_pss,C_pss,D_pss] = ssdata(net.a_controller_local{1}.pss.get_sys);

%1axis, num_x = 7
%size of A(10,10)
A = [0 omega0 zeros(1,8);
     0 -D/M   zeros(1,8);
     zeros(4,10);
     0  D_pss*(Ka/Te)  0 0 0 0 -1/Te C_pss.*(Ka/Te);
     [zeros(3,1),B_pss,zeros(3,5),A_pss]];

%size of L(10,9)
L = [zeros(1,9);
     -1/M zeros(1,7) 1/M;
     0 -1/Tdo 0 0 0 0 1/Tdo 0 0;
     0  0 -1/Tqo 0 0 0 0 0 0;
     0 0 0 1/TTqo 0 0 0 0 0;
     0 0 0 0 1/TTdo 0 0 0 0;
     0 0 0 0 0 -Ka/Te 0 -Ka/Te 0;
     zeros(3,9);
     ];

%size of gamma(7,10)
gamma = [eye(7,7),zeros(7,3)];


%too lazy to differentiate v with w
%Define syms
syms s_delta s_omega s_Eq s_Ed s_psiq s_psid s_Vap s_Vpss s_Pmech
%Define Vd, Vq, Id, Iq
Vq = Vabs*cos(s_delta-Vangle);
Vd = Vabs*sin(s_delta-Vangle);

for_Id = (Xdpp-Xls)*s_Eq/(Xdp-Xls) + (Xdp-Xdpp)*s_psid/(Xdp-Xls);
for_Iq = (Xqpp-Xls)*s_Ed/(Xqp-Xls) - (Xqp-Xqpp)*s_psiq/(Xqp-Xls);
Iq = (Vd - for_Iq)/Xqpp;
Id = (for_Id - Vq)/Xdpp;

P = Vq*Iq+Vd*Id;
Efd = s_Eq + (Xd-Xdp)*(Id- ((Xdp-Xdpp)/(Xdp-Xls)^2)*(s_psid+(Xdp-Xls)*Id-s_Eq));
Efq = s_Ed - (Xq-Xqp)*(Iq- ((Xqp-Xqpp)/(Xqp-Xls)^2)*(s_psiq+(Xqp-Xls)*Iq+s_Ed));
psiq_ = -s_psiq-s_Ed-(Xqp-Xls)*Iq;
psid_ = -s_psid+s_Eq-(Xdp-Xls)*Id;
Vfd = s_Vap;

Environ = [ P; Efd; Efq; psiq_; psid_; Vabs; Vfd; s_Vpss; s_Pmech];
J = jacobian(Environ, [s_delta, s_omega, s_Eq, s_Ed, s_psiq, s_psid, s_Vap]);
environ = eval(subs(J,[s_delta,s_Eq,s_Ed,s_psiq,s_psid],[delta,Eq,Ed,psiq,psid]));
A_ = (L*environ*gamma);
Ahat = A+A_;
Bw = L*environ;




