net = network_IEEE68bus;

%Find generator index
bus_num    = numel(net.a_bus); 
component_idx = zeros(1,bus_num);

for idx = 1:bus_num
    switch class(net.a_bus{idx}.component) %can return class names of buses
        case 'generator_1axis'
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

for i = 1:length(generator_idx)
    retrofit_hinf = controller_retrofit_hinf(net,generator_idx(i));
    net.add_controller_local(retrofit_hinf);
end

%Controller AGC
con = controller_broadcast_PI_AGC(net,1:2,1:2,-10,-500);
net.add_controller_global(con);

%fault
option.fault = {{[0 0.01], 1}, {[15 15.01], 1}, {[10,10.01],2}};

time = [0,10,20,60];
out = net.simulate(time,option);

for i = 1:length(generator_idx)
    plot(out.t,out.X{i}(:,2));
    hold on;
end