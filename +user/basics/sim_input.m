% ネットワークの定義
net = network_sample3bus();
%条件設定
time = [0,10,20,60];
%input to bus number 3, load
u_idx = 3;
%
u = [0, 0.05, 0.1, 0.1;...
     0,    0,   0,   0];

%解析実行
out = net.simulate(time,u, u_idx);
plot(out.t,[out.X{1}(:,2),out.X{2}(:,2)])
