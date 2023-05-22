net = network_IEEE68bus;

branch_num = numel(net.a_branch); %68 buses
bus_num    = numel(net.a_bus); %83 cells

adjacency_matrix = zeros(bus_num,bus_num);

for idx = 1:branch_num
    from = net.a_branch{idx}.from;
    to   = net.a_branch{idx}.to;
    adjacency_matrix(from,to)=1;
    adjacency_matrix(to,from)=1;
end

component_idx = zeros(1,bus_num);
%Depends on class names of components, assign numbers
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

h = plot(graph(adjacency_matrix));
highlight(h,load_idx,Marker="^", MarkerSize=5)
highlight(h,generator_idx,Marker='o',MarkerSize=5)
highlight(h,empty_idx,Marker = ".",MarkerSize=5)
