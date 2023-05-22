function init_park
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
component1 = generator_PARK(omega0,machinery(1,:));
%component1 = generator_1axis(omega0,mac_data);
net.a_bus{1}.set_component(component1);

Xd = 1.220; Xd_prime = 0.667; Xq = 0.667; Tdo = 8.97; M = 12; D = 10;
mac_data2 = table(Xd,Xd_prime,Xq,Tdo,M,D);
component2 = generator_PARK(omega0, machinery(2,:));
%component2 = generator_1axis(omega0, mac_data2);
net.a_bus{2}.set_component(component2);

component3 = load_impedance();
net.a_bus{3}.set_component(component3);

%AGC
con = controller_broadcast_PI_AGC(net,1:2,1:2,-10,-500);
net.add_controller_global(con);

net.initialize;

option          = struct();

time = [0,100];
option.fault = {{[10 10.08], 1}};
out = net.simulate(time, option);
plot(out.t,[out.X{1}(:,2)]);
end