net = network_IEEE68bus;
bus_num = numel(net.a_bus);

%arrayfun(@(x)func(x), A) means apply func(x) to every element of A
component_list = arrayfun(@(idx) {['bus',num2str(idx)] , class(net.a_bus{idx}.component)},(1:bus_num)','UniformOutput',false);
cell2table(vertcat(component_list{:}) ,"VariableNames",["idx" "component"])



