function axis1_hinf
%Create an empty power network first
net = power_network();

%BRANCH
%branch_pi(from,to,[xreal, ximag],b) b is ground capacitance
branch12 = branch_pi(1,2,[0.010, 0.085],0);
net.add_branch(branch12);

branch23 = branch_pi(2,3,[0.017,0.092],0);
net.add_branch(branch23);

%BUS
%bus_PQ, bus_PV, bus_slack(Vabs,Vangle,shunt)
shunt = [0,0];
bus_1 = bus_slack(2,0,shunt);
net.add_bus(bus_1);
bus_2 = bus_PV(0.5,2,shunt);
net.add_bus(bus_2);
bus_3 = bus_PQ(-3,0,shunt);
net.add_bus(bus_3);

%COMPONENT
%component_type_name(omega0, mac_data)
omega0 = 60*2*pi;
Xd = 1.569; Xd_prime = 0.963; Xq = 0.963; Tdo = 5.14; M = 100; D = 10;
mac_data = table(Xd,Xd_prime,Xq,Tdo,M,D);
machinery = readtable('+user/3bus_PARK/machinery.xlsx');
component1 = generator_1axis(omega0,machinery(1,:));
%component1 = generator_PARK(omega0,mac_data);
net.a_bus{1}.set_component(component1);

Xd = 1.220; Xd_prime = 0.667; Xq = 0.667; Tdo = 8.97; M = 12; D = 10;
mac_data2 = table(Xd,Xd_prime,Xq,Tdo,M,D);
%component2 = generator_PARK( omega0, mac_data2);
component2 = generator_1axis(omega0, machinery(2,:));
net.a_bus{2}.set_component(component2);

component3 = load_impedance();
net.a_bus{3}.set_component(component3);

%Set PSS
No_bus = 1; Kpss = 250; Tpss = 10; TL1p = 0.07; TL1 = 0.02; TL2p = 0.07; TL2 = 0.02;
pss_data = table(No_bus,Kpss,Tpss,TL1p,TL1,TL2p,TL2);
net.a_bus{1}.component.set_pss(pss(pss_data));
net.a_bus{2}.component.set_pss(pss(pss_data));

%Set AVR
No_bus = 1; Ka = 20; Te = 0.05;
avr_data = table(No_bus,Ka,Te);
net.a_bus{1}.component.set_avr(avr_sadamoto2019(avr_data));
net.a_bus{2}.component.set_avr(avr_sadamoto2019(avr_data));

%Initialize power flow calculation to store new properties in power_network
net.initialize;
%get admittance matrix
net.get_admittance_matrix;
full(net.get_admittance_matrix);

%Controller AGC

con = controller_broadcast_PI_AGC(net,1:2,1:2,-10,-500);
net.add_controller_global(con);

%Controller retrofit LQR

idx = 1;
LQR = controller_retrofit_hinf(net,idx);
net.add_controller_local(LQR);
idx = 2;
LQR2 = controller_retrofit_hinf(net,idx);
net.add_controller_local(LQR2);

%{
idx = 1;
retrofit_hinf1 = controller_retrofit_hinf(net,idx);
net.add_controller_local(retrofit_hinf1);
idx = 2;
retrofit_hinf2 = controller_retrofit_hinf(net,idx);
net.add_controller_local(retrofit_hinf2);
%}

%simulate

%initial value change
option          = struct();
option.x0_sys  = net.x_equilibrium;
%option.x0_sys(2)= option.x0_sys(2) + 0.001;


%load change, input
%{
u_idx = 3;
u     = [0, 0.05, 0.1, 0.1;...
         0,    0,   0, 0];
%}


%fault
time = [0,100];
option.fault = {{[10 10.01], 1}};
out = net.simulate(time,option);
plot(out.t,out.U{1}(:,1));

end