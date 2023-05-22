net = network_3bus_1axis();
net.simulate

function net = network_3bus_1axis()

omega0 = 60*2*pi;
net = power_network();

bus_array = [1 1.00 0 0.60 0 0 0 0 0 2;
             2 1.00 0 5.45 0 1.0 1.0 0 0 3;
             3 1.00 0 1.00 2 0 0 0 0 1];
bus = array2table(bus_array, 'VariableNames', ...
    {'No', 'V_abs', 'V_angle', 'P_gen', 'Q_gen', 'P_load', 'Q_load', 'G_shunt', 'B_shunt', 'type'} ...
    );

branch_array = [...
            1 2 0   1/12.56041  0   1   0;...
            3 2 0   1/13.65107  0   1   0;];
branch = array2table(branch_array, 'VariableNames', ...
                {'bus_from', 'bus_to', 'x_real', 'x_imag', 'y', 'tap', 'phase'}...
                );

mac_array = [1, 1, 0.1,    0.031,  0.069, 10.2, 84,   4;
                   2, 3, 0.295,  0.0697, 0.282, 6.56, 60.4, 9.75];
mac_data = array2table(mac_array, 'VariableNames', ...
            {'No_machine', 'No_bus', 'Xd', 'Xd_prime', 'Xq', 'T', 'M', 'D'} ...
            );

exc_array = [1 0 0.05;
                    3 0 0.05];
exc_data = array2table(exc_array, 'VariableNames', ...
            {'No_bus', 'Ka', 'Te'} ...
            );

pss_array = [1 0 10 0.05 0.015 0.08 0.01;
             3 0 10 0.05 0.015 0.08 0.01];
pss_data = array2table(pss_array, 'VariableNames', ...
            {'No_bus', 'Kpss', 'Tpss', 'TL1p', 'TL1', 'TL2p', 'TL2'} ...
            );

for i = 1:size(bus, 1)
    shunt = bus{i, {'G_shunt', 'B_shunt'}};
    switch bus{i, 'type'}
        case 1
            V_abs = bus{i, 'V_abs'};
            V_angle = bus{i, 'V_angle'};
            b = bus_slack(V_abs, V_angle, shunt);
            b.set_component(get_generator(i, mac_data, exc_data, pss_data, omega0));

        case 2
            V_abs = bus{i, 'V_abs'};
            P = bus{i, 'P_gen'};
            b = bus_PV(P, V_abs, shunt);
            b.set_component(get_generator(i, mac_data, exc_data, pss_data, omega0));

        case 3
            P = bus{i, 'P_load'};
            Q = bus{i, 'Q_load'};
            b = bus_PQ(-P, -Q, shunt);
            if P~=0 || Q~=0
                load = load_impedance();
                b.set_component(load);
            end

    end
    net.add_bus(b);
end

for i = 1:size(branch, 1)
    if branch{i, 'tap'} == 0
        br = branch_pi(branch{i, 'bus_from'}, branch{i, 'bus_to'},...
            branch{i, {'x_real', 'x_imag'}}, branch{i, 'y'});
    else
        br = branch_pi_transformer(branch{i, 'bus_from'}, branch{i, 'bus_to'},...
            branch{i, {'x_real', 'x_imag'}}, branch{i, 'y'},...
            branch{i, 'tap'}, branch{i, 'phase'});
    end
    net.add_branch(br);
end

net.initialize();
end


function g = get_generator(i, mac_data, exc_data, pss_data, omega0)
idx = mac_data{:, 'No_bus'} == i;
if sum(idx) ~= 0
    g = generator_1axis(omega0, mac_data(idx, :));
    exc = exc_data(exc_data{:, 'No_bus'}==i, :);
    g.set_avr(avr_sadamoto2019(exc));
    p = pss_data(pss_data{:, 'No_bus'}==i, :);
    g.set_pss(pss(p));
end
end
