function park_hinf
net = park_IEEE68bus;

%Find generator index
bus_num    = numel(net.a_bus); 
component_idx = zeros(1,bus_num);

for idx = 1:bus_num
    switch class(net.a_bus{idx}.component) %can return class names of buses
        case 'generator_PARK'
            component_idx(idx) = 1;
        case 'load_impedance'
            component_idx(idx) = 2;
        case 'component_empty'
            component_idx(idx) = 3;
    end
end

generator_idx = find(component_idx==1);
load_idx = find(component_idx==2);
empty_idx = find(component_idx==3);
n = length(generator_idx);

%AGC
con = controller_broadcast_PI_AGC(net,1:n,1:n,-10,-500);
net.add_controller_global(con);


for i = 1:n
    LQR = controller_retrofit_hinf_PARK(net,i);
    net.add_controller_local(LQR);
end


%initial value change
option          = struct();
option.x0_sys  = net.x_equilibrium;
option.x0_sys(2)= option.x0_sys(2) + 0.001;
%load change, input
u_idx = 25;
u     = [0, 0.05, 0.1, 0.1;...
         0,    0,   0, 0];

%fault
option.fault = {{[0 0.01], 8},{[20 20.01], 8},{[40 40.01], 8}};

time = [0,10,30,60];
out = net.simulate(time,u,u_idx, option);
plot(out.t,[out.X{8}(:,2)]);
end