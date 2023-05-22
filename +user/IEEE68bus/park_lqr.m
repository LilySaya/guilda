
net = park_IEEE68bus;


k = 1e-5;
 Q  = [k, 1e3, k, k,   k, k, k,k,k,k];           
 R  = 1e-9;
 
%AGC
con = controller_broadcast_PI_AGC(net,1:16,1:16,-10,-500);
net.add_controller_global(con);


for i = 1:16
    LQR = controller_retrofit_LQR_PARK(net,i,diag(Q),diag(R));
    net.add_controller_local(LQR);
end
 

%initial value change
option          = struct();
option.x0_sys  = net.x_equilibrium;
%option.x0_sys(2)= option.x0_sys(2) + 0.001;
%load change, input
%{
u_idx = 25;
u     = [0, 0.05, 0.1, 0.1;...
         0,    0,   0, 0];
%}
%fault
%{
option.fault = {{[0 0.1], 1}};
out3 = net.simulate([0 60],option);
%plot(out.t,out.X{11}(:,2));
%}

out3 = net.simulate([0,60],...
                    'fault',     {{[0,0.1],1}},...
                    'OutputFcn',  'omega', ...
                    'tools',     false, ...
                    'do_report', false);

calculate_l2norm(out3);
calculate_settlingtime(out3);

function T = calculate_settlingtime(out)
    T = zeros(1,16);
    Peak = zeros(1,16);
    for i = 1:16
        x_i = out.X{i}(:,2);
        t_i = out.t;
        info = stepinfo(x_i,t_i,0);
        Peak(i) = info.Peak;
    end
    P = max(Peak);
    clear info;


    disp(' ')
    disp(['最大値ピーク値は',num2str(P)])
    %Pを自分で設定し固定することもできる
    %P = 0.002;
    error = 0.02;
    disp(['過渡時間:周波数偏差の絶対値が',num2str(P*error),'以下になる時間'])
    fprintf('\n各発電機の周波数偏差の過渡時間(大きい順) \n\n')
    for i = 1:16
        x_i = out.X{i}(:,2);
        t_i = out.t;
        %過渡時間のエラーwhich depends on each peakを統一したいのでmax Peakを使う
        info =  stepinfo(x_i,t_i,0,'SettlingTimeThreshold',P*error/Peak(i)); %定常値0
        T(i) = info.TransientTime;
        clear info;
    end
    [~, idx_sort] = sort(T,'descend');

    for i=idx_sort
        disp(['同期発電機',num2str(i,'%.2d'),'の過渡時間: ',num2str(T(i))])      
    end
    disp(' ')
    disp(['全発電機の過渡時間の和：',num2str(sum(T))])
end


function L2 = calculate_l2norm(out)
    L2 = zeros(1,16);                       %データ格納用の行列を定義
    fprintf('\n各発電機の周波数偏差のL２ノルム(大きい順) \n\n')
    for i = 1:16                            %IEEE68busは母線1~16に同期発電機が付加されている．
        omega_i  = out.X{i}(:,2);     %母線iの周波数偏差の応答データを抽出
        omega2_i = omega_i.^2;              %応答データの各要素を２乗する
        integer  = trapz(out.t, omega2_i);  %２乗したデータを数値積分する
        L2_i     = sqrt(integer);           %積分値の平方根を取りL２ノルムを計算
        L2(i)    = L2_i;                    %得られたL2ノルムを配列に格納する
    end
    [~, idx_sort] = sort(L2,'descend');     %L2ノルムの大きい順にソートした機器の番号を得る
    
    for i = idx_sort
        disp(['同期発電機',num2str(i,'%.2d'),'のl2ノルム: ',num2str(L2(i))])
    end
    disp(' ')
    disp(['全発電機のL２ノルムの和：',num2str(sum(L2))])
end
