function L2 = L2_norm(out)
    L2 = zeros(1,16);                       %データ格納用の行列を定義
    fprintf('\n各発電機の周波数偏差のL２ノルム(大きい順) \n\n')
    for i = 1:16                            %IEEE68busは母線1~16に同期発電機が付加されている．
        omega_i  = out.X{i}{:,'omega'};     %母線iの周波数偏差の応答データを抽出
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