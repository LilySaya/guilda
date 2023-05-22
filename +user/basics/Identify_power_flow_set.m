net = network_IEEE68bus;
bus_num = numel(net.a_bus);
%prepare an empty cell with number of buses
bus_idx = cell(bus_num,1);

Vabs  = nan(bus_num,1);
Vangle= nan(bus_num,1);
P     = nan(bus_num,1);
Q     = nan(bus_num,1);

for idx = 1:bus_num
    switch class(net.a_bus{idx})
        case 'bus_PV'
            P(idx)     = net.a_bus{idx}.P;
            Vabs(idx)  = net.a_bus{idx}.Vabs;
        case 'bus_PQ'
            P(idx)     = net.a_bus{idx}.P;
            Q(idx)     = net.a_bus{idx}.Q;
        case 'bus_slack'
            Vabs(idx)  = net.a_bus{idx}.Vabs;
            Vangle(idx)= net.a_bus{idx}.Vangle;
    end
end
idx = (1:bus_num)';
powerflow_set = table(idx,Vabs,Vangle,P,Q)



