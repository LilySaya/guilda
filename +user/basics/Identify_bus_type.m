%define network
net = network_IEEE68bus;
%get number of buses
bus_num = numel(net.a_bus);
%Prepare list of 0 with number of buses
bus_idx = zeros(1,bus_num);

%Depends on class names of buses, assign numbers
for idx = 1:bus_num
    switch class(net.a_bus{idx}) %can return class names of buses
        case 'bus_PV'
            bus_idx(idx) = 1;
        case 'bus_PQ'
            bus_idx(idx) = 2;
        case 'bus_slack'
            bus_idx(idx) = 3;
    end
end

%Find corresponding numbers and print out where they are
PV_bus_idx = find(bus_idx==1)
PQ_bus_idx = find(bus_idx==2)
slack_bus_idx = find(bus_idx==3)
