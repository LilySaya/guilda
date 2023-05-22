% ネットワークの定義
net = network_IEEE68bus;

time = [0,10,20,60];

%Initial_value_response
bus_idx   = 3; %3番目の母線
state_idx = 2; %2番目の状態
%母線1から母線(bus_idx-1)番目までの状態の個数を足し合わせる
idx = 0;
for i = 1:bus_idx-1
    idx = idx+net.a_bus{i}.component.get_nx;
end
%状態の2番目のインデックスを知りたいため
idx = idx+state_idx;
%IEEE68busモデルの場合は、idx = 7+7+2= 16となります。

option = struct();
option.x0_sys = net.x_equilibrium;
option.x0_sys(idx)= option.x0_sys(idx) + 0.01;

%Fault
option.fault = {{[0 0.07], 1}, {[15 15.05], 1}, {[10,10.01],[2,5:7]}};


%Input
u_idx = 20;

u = [0, 0.05, 0.1, 0.1;...
     0,    0,   0,   0];

%SimulationResult
option.tools = true;
%解析実行
out = net.simulate(time,u, u_idx,option);
%plot(out.t,[out.X{1}(:,2),out.X{2}(:,2)]);

%SimulationResult UI
out.UIplot