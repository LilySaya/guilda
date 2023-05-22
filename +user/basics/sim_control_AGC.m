time = [0,10,20,60];
u_idx = 3;
u     = [0, 0.05, 0.1, 0.1;...
         0,    0,   0,   0];


%AGCコントローラを定義
con = controller_broadcast_PI_AGC(net,1:2,1:2,-10,-50);

%電力系統にcontrollerクラスを代入
net.add_controller_global(con);

out1 = net.simulate(time,u, u_idx);

%データ抽出
sampling_time = out1.t;
omega1 = out1.X{1}(:,2);
omega2 = out1.X{2}(:,2);

%プロット
figure;
hold on;
plot(sampling_time, omega2,'LineWidth',2)
plot(sampling_time, omega1,'LineWidth',2)
xlabel('時刻(s)','FontSize',15);
ylabel('周波数偏差','FontSize',15);
legend({'機器2の周波数偏差','機器1の周波数偏差'})
title('各同期発電機の周波数偏差','FontSize',20)
hold off

