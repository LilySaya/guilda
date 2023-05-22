% ネットワークの定義
net = network_IEEE68bus();

% シミュレーションのためのオプションを定義・決定
option          = struct();
option.x0_sys  = net.x_equilibrium;
option.x0_sys(16)= option.x0_sys(16) + 0.01;

% シミュレーションの実行
%out = net.simulate([0 20], option);

%without using struct
x0    = net.x_equilibrium;
x0(16)= x0(16) + 0.01;
%_out = net.simulate([0 20],'x0_sys',x0);

%using get_nx
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

option.x0_sys(idx)= option.x0_sys(idx) + 0.1;
out = net.simulate([0 20],option)           

