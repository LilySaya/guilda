classdef generator_PARK < component % 状態・パラメーターはqを先においている
    %edited
    %VとIの関係をedit
    properties(Access = private)
        parameter_vec
        system_matrix
        x_st
        V_st
        I_st
    end
    
    properties(SetAccess = private)
         parameter
         x_equilibrium
         V_equilibrium
         I_equilibrium
        avr
        pss
        governor
        alpha_st
        omega0
    end
    
    methods
        function obj = generator_PARK(omega, parameter)
            obj.omega0 = omega;
            if isstruct(parameter)
                parameter = struct2table(parameter);
            end
            % PARK用のパラメータ名に変更
            obj.parameter = parameter(:, {'Xq', 'Xq_prime', 'Xq_pp','Xd', 'Xd_prime', 'Xd_pp','X_ls','Tdo', 'Tqo', 'TTdo','TTqo','M', 'D'});   % ソートしてるだけ
            obj.parameter_vec = obj.parameter.Variables;
            obj.avr = avr();
            obj.governor = governor();
            obj.pss = pss();
            obj.system_matrix = struct();
        end
        
        function name_tag = get_x_name(obj)
            %Added psiq, psid
            gen_state = {'delta','omega','Eq','Ed','psiq','psid'};
            avr_state = obj.avr.get_state_name;
            pss_state = obj.pss.get_state_name;
            governor_state = obj.governor.get_state_name;
            name_tag = horzcat(gen_state,avr_state,pss_state,governor_state);
        end

        function u_name = get_u_name(obj)
            u_name = {'Vfd','Pm'};
        end
        
        function out = get_nx(obj)
            % 右辺の第１項を4から6に変更
            out = 6 + obj.avr.get_nx() + obj.pss.get_nx() + obj.governor.get_nx();
        end
        
        function nu = get_nu(obj)
            nu = 2;
        end
        
        function [dx, con] = get_dx_constraint(obj, t, x, V, I, u)
            % このx,V,Iは平衡状態ではなく、その時刻における値
            omega0 = obj.omega0;
            % パラメータを追加
            Xq = obj.parameter_vec(1);
            Xqp = obj.parameter_vec(2);
            Xqpp = obj.parameter_vec(3);
            Xd = obj.parameter_vec(4);
            Xdp = obj.parameter_vec(5);
            Xdpp = obj.parameter_vec(6);
            Xls = obj.parameter_vec(7);
            Tdo = obj.parameter_vec(8);
            Tqo = obj.parameter_vec(9);
            TTdo = obj.parameter_vec(10);
            TTqo = obj.parameter_vec(11);
            M = obj.parameter_vec(12);
            d = obj.parameter_vec(13);
            % nxの値を4から6に変更
            nx = 6;
            nx_avr = obj.avr.get_nx();
            nx_pss = obj.pss.get_nx();
            nx_gov = obj.governor.get_nx();
            
            x_gen = x(1:nx);
            x_avr = x(nx+(1:nx_avr));
            x_pss = x(nx+nx_avr+(1:nx_pss));
            x_gov = x(nx+nx_avr+nx_pss+(1:nx_gov));
            
            Vabs = norm(V);
            Vangle = atan2(V(2), V(1));
            
            delta = x_gen(1);
            omega = x_gen(2);
            Eq = x_gen(3);
            Ed = x_gen(4);
            psiq = x_gen(5);
            psid = x_gen(6);
            
            % Vd, Vqを定義
            Vq = V(1)*cos(delta)+V(2)*sin(delta); %Vabs*cos(delta-Vang)
            Vd = V(1)*sin(delta)-V(2)*cos(delta); %Vabs*sin(delta-Vang)
            
            %omega を無視せずにId,Iq を定義
            for_Id = (1+omega)*((Xdpp-Xls)*Eq/(Xdp-Xls) + (Xdp-Xdpp)*psid/(Xdp-Xls));
            for_Iq = (1+omega)*((Xqpp-Xls)*Ed/(Xqp-Xls) - (Xqp-Xqpp)*psiq/(Xqp-Xls));
            Iq = (Vd - for_Iq)/(Xqpp);
            Id = (for_Id - Vq)/(Xdpp);
            %Iq = (Vd - for_Iq)/((1+omega)*Xqpp);
            %Id = (for_Id - Vq)/((1+omega)*Xdpp);
            
        
            
            % Id, Iqを定義のために
            %{
            for_Id = (Xdpp-Xls)*Eq/(Xdp-Xls) + (Xdp-Xdpp)*psid/(Xdp-Xls);
            for_Iq = (Xqpp-Xls)*Ed/(Xqp-Xls) - (Xqp-Xqpp)*psiq/(Xqp-Xls);
            Iq = (Vd - for_Iq)/Xqpp;
            Id = (for_Id - Vq)/Xdpp;
            %}

            % Efdの修正とEfqの追加
            Efd = Eq + (Xd-Xdp)*(Id- ((Xdp-Xdpp)/(Xdp-Xls)^2)*(psid+(Xdp-Xls)*Id-Eq));
            Efq = Ed - (Xq-Xqp)*(Iq- ((Xqp-Xqpp)/(Xqp-Xls)^2)*(psiq+(Xqp-Xls)*Iq+Ed));
            
            [dx_pss, v] = obj.pss.get_u(x_pss, omega);
            [dx_avr, Vfd] = obj.avr.get_Vfd(x_avr, Vabs, Efd, u(1)-v);
            [dx_gov, P] = obj.governor.get_P(x_gov, omega, u(2));
            
            % dEをdEqに，dEdの追加
            dEq = (-Efd + Vfd)/Tdo;
            dEd = (-Efq)/Tqo;
            ddelta = obj.omega0 * omega; 
            % PはPmechを指す
            domega = (P - d*omega - Vq*Iq - Vd*Id)/M;

            %dpsid, dpsiqを追加
            psiq_ = -psiq-Ed-(Xqp-Xls)*Iq;
            psid_ = -psid+Eq-(Xdp-Xls)*Id;
            dpsiq = psiq_/TTqo;
            dpsid = psid_/TTdo;

            % |I|cosI, |I|sinIを逆算
            Ir = Id*sin(delta)+Iq*cos(delta);
            Ii = -Id*cos(delta)+Iq*sin(delta);
            
            con = I - [Ir; Ii];
            % ここで，dEqとdEdの順序を逆にする必要があるかも
            dx = [ddelta; domega; dEq; dEd; dpsiq; dpsid; dx_avr; dx_pss; dx_gov];
        end %ok
        
        function [dx, con] = get_dx_constraint_linear(obj, t, x, V, I, u)
            A  = obj.system_matrix.A;
            B  = obj.system_matrix.B;
            C  = obj.system_matrix.C;
            D  = obj.system_matrix.D;
            BV = obj.system_matrix.BV;
            DV = obj.system_matrix.DV;
            BI = obj.system_matrix.BI;
            DI = obj.system_matrix.DI;
            dx = A*(x-obj.x_st) + B*u + BV*(V-obj.V_st) + BI*(I-obj.I_st);
            con = C*(x-obj.x_st) + D*u + DV*(V-obj.V_st) + DI*(I-obj.I_st);
        end

        function [A, B, C, D, BV, DV, BI, DI, R, S] = get_linear_matrix(obj, x_st, Vst)
            if nargin < 2 || (isempty(x_st) && isempty(Vst))
                A  = obj.system_matrix.A;
                B  = obj.system_matrix.B;
                C  = obj.system_matrix.C;
                D  = obj.system_matrix.D;
                BV = obj.system_matrix.BV;
                DV = obj.system_matrix.DV;
                BI = obj.system_matrix.BI;
                DI = obj.system_matrix.DI;
                R = obj.system_matrix.R;
                S = obj.system_matrix.S;
                return;
            end
            if nargin < 2 || isempty(x_st)
                x_st = obj.x_st;
            end
            if nargin < 3 || isempty(Vst)
                Vst = obj.V_st;
            end
            omega0 = obj.omega0;
            % パラメータを追加
            Xq = obj.parameter_vec(1);
            Xqp = obj.parameter_vec(2);
            Xqpp = obj.parameter_vec(3);
            Xd = obj.parameter_vec(4);
            Xdp = obj.parameter_vec(5);
            Xdpp = obj.parameter_vec(6);
            Xls = obj.parameter_vec(7);
            Tdo = obj.parameter_vec(8);
            Tqo = obj.parameter_vec(9);
            TTdo = obj.parameter_vec(10);
            TTqo = obj.parameter_vec(11);
            M = obj.parameter_vec(12);
            d = obj.parameter_vec(13);
            
            % 表記を変更：Vabscos -> Vq, Vabssin -> Vd

            % 動揺方程式の状態変数を４つに変更
            % x1 = delta
            % x2 = omega
            % x3 = Eq
            % x4 = Ed
            % x5 = psiq
            % x6 = psid
            A_swing = [0 omega0 0 0 0 0;
                0 -d/M 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0
                0 0 0 0 0 0
                0 0 0 0 0 0];
            % u1 = Pmech
            % u2 = Vfd
            % u3 = Pout
            % u4 = Efd
            % u5 = Efq
            % u6 = psiq_
            % u7 = psid_
            B_swing = [0, 0, 0, 0, 0, 0, 0;
                1/M, 0, -1/M, 0, 0, 0, 0;
                0, 1/Tdo, 0, -1/Tdo, 0, 0, 0;
                0, 0, 0, 0, -1/Tqo, 0, 0;
                0, 0, 0, 0, 0, 1/TTqo, 0;
                0, 0, 0, 0, 0, 0, 1/TTdo];
            % y = [delta, omega, Eq, Ed, psiq, psid]
            C_swing = eye(6);
            sys_swing = ss(A_swing, B_swing, C_swing, 0);
            OutputGroup = struct();
            OutputGroup.delta = 1;
            OutputGroup.omega = 2;
            OutputGroup.Eq = 3;
            OutputGroup.Ed = 4;
            OutputGroup.psiq = 5;
            OutputGroup.psid = 6;
            sys_swing.OutputGroup = OutputGroup;
            InputGroup = struct();
            InputGroup.Pmech = 1;
            InputGroup.Vfd = 2;
            InputGroup.Pout = 3;
            InputGroup.Efd_swing = 4;
            InputGroup.Efq_swing = 5;
            InputGroup.psiq_ = 6;
            InputGroup.psid_ = 7;
            sys_swing.InputGroup = InputGroup;
            
            % テイラー展開を用いるので，ここからは平衡状態
            delta = x_st(1);
            Eq = x_st(3);
            Ed = x_st(4);
            psiq = x_st(5);
            psid = x_st(6);
            
            syms s_delta s_Eq s_Ed s_psiq s_psid s_Vreal s_Vimag
            Vq = s_Vreal*cos(s_delta)+s_Vimag*sin(s_delta);%Vabs*cos(delta-Vang)
            Vd = s_Vreal*sin(s_delta)-s_Vimag*cos(s_delta);%Vabs*sin(delta-Vang)

            for_Id = (Xdpp-Xls)*s_Eq/(Xdp-Xls) + (Xdp-Xdpp)*s_psid/(Xdp-Xls);
            for_Iq = -(Xqpp-Xls)*s_Ed/(Xqp-Xls) + (Xqp-Xqpp)*s_psiq/(Xqp-Xls);
            Id = (for_Id-Vq)/Xdpp;
            Iq = (for_Iq+Vd)/Xqpp;
            
            % |I|cosI, |I|sinIを逆算
            Ir = Id*sin(s_delta)+Iq*cos(s_delta);
            Ii = -Id*cos(s_delta)+Iq*sin(s_delta);
            I = [Ir;Ii]; %1%2

            P = [s_Vreal, s_Vimag]*I; %3
            Efd = s_Eq + (Xd-Xdp)*(Id- ((Xdp-Xdpp)/(Xdp-Xls)^2)*(s_psid+(Xdp-Xls)*Id-s_Eq));%4
            Efq = s_Ed - (Xq-Xqp)*(Iq- ((Xqp-Xqpp)/(Xqp-Xls)^2)*(s_psiq+(Xqp-Xls)*Iq+s_Ed));%5
            psiq_ = -s_psiq-s_Ed-(Xqp-Xls)*Iq; %6
            psid_ = -s_psid+s_Eq-(Xdp-Xls)*Id; %7

            J = jacobian([P, Efd, Efq, psiq_, psid_, I'], ...
                      [s_delta, s_Eq, s_Ed, s_psiq,s_psid, s_Vreal, s_Vimag]);
            sys_fb = ss(eval(subs(J,[s_delta, s_Eq, s_Ed, s_psiq, s_psid, s_Vreal, s_Vimag],[delta, Eq, Ed, psiq, psid, Vst(1), Vst(2)])));
            InputGroup = struct();
            InputGroup.delta = 1;
            InputGroup.Eq = 2;
            InputGroup.Ed = 3;
            InputGroup.psiq = 4;
            InputGroup.psid = 5;
            InputGroup.V = 6:7;
            sys_fb.InputGroup = InputGroup;
            OutputGroup = struct();
            OutputGroup.P = 1;
            OutputGroup.Efd = 2;
            OutputGroup.Efq = 3;
            OutputGroup.psiq_ = 4;
            OutputGroup.psid_ = 5;
            OutputGroup.I = 6:7;
            sys_fb.OutputGroup = OutputGroup; %ok
            
            % これも平衡点
            Vabs = norm(Vst);
            
            sys_V = ss([eye(2); Vst'/Vabs]);
            sys_V.InputGroup.Vin = 1:2;
            OutputGroup = struct();
            OutputGroup.V = 1:2;
            OutputGroup.Vabs = 3;
            sys_V.OutputGroup = OutputGroup; %ok
            
            sys_avr = obj.avr.get_sys(); %avrのシステムを取得
            sys_pss = obj.pss.get_sys(); %pssのシステムを取得
            sys_gov = obj.governor.get_sys(); %governorのシステムを取得
            G = blkdiag(sys_swing, sys_fb, sys_V, sys_avr, -sys_pss, sys_gov);
            ig = G.InputGroup;
            og = G.OutputGroup;
            % Efqを追加、Eq,Edに変更、ig.Vinとog.Vをつなげた
            feedin = [ig.Pout, ig.Efd, ig.Efd_swing, ig.Efq_swing, ig.delta, ig.Eq, ig.Ed, ig.V, ig.Vabs, ig.Vfd, ig.u_avr, ig.omega, ig.omega_governor, ig.Pmech];
            feedout = [og.P, og.Efd, og.Efd, og.Efq, og.delta, og.Eq, og.Ed, og.V, og.Vabs, og.Vfd, og.v_pss, og.omega, og.omega, og.Pmech];
            I = ss(eye(numel(feedin)));
            
            ret = feedback(G, I, feedin, feedout, 1);

            % ここから先はよくわかっていない
            ret_u = ret('I', {'u_avr',  'u_governor'});
            ret_V = ret('I', 'Vin');
            A = ret.a;
            B = ret_u.b;
            C = ret_u.c;
            D = ret_u.d;
            BV = ret_V.b;
            DV = ret_V.d;
            BI = zeros(size(A, 1), 2);
            DI = -eye(2);
            R = BV;
            S = zeros(1, size(A, 1));
            S(2) = 1;
        end
        
        function set_avr(obj, avr)
            if isa(avr, 'avr')
                obj.avr = avr;
            else
               error(''); 
            end
        end
        
        function set_pss(obj, pss)
            if isa(pss, 'pss')
                obj.pss = pss;
            else
                error('');
            end
        end

        function set_governor(obj, governor)
            if isa(governor, 'governor')
                obj.governor = governor;
            else
                error('');
            end
        end
        
        function initialize_net(obj)
            if ~isempty(obj.net)
                obj.net.initialize(false);
            end
        end
        
        function set_linear_matrix(obj, varargin)
            if isempty(obj.omega0)
                return
            end
            mat = struct();
            [mat.A, mat.B, mat.C, mat.D, mat.BV, mat.DV, mat.BI, mat.DI, mat.R, mat.S] = obj.get_linear_matrix(varargin{:});
            obj.system_matrix = mat;
        end
        
        % 潮流計算結果から逆算して平衡点を算出
        function x_st = set_equilibrium(obj, V, I)
            Vangle = angle(V);
            Vabs = abs(V);
            Iangle = angle(I);
            Iabs = abs(I);
            Pow = conj(I)*V;
            P = real(Pow);
            Q = imag(Pow);
            Xd = obj.parameter{:, 'Xd'};
            Xdp = obj.parameter{:, 'Xd_prime'};
            Xq = obj.parameter{:, 'Xq'};
            Xqp = obj.parameter{:, 'Xq_prime'};
            Xls = obj.parameter{:,'X_ls'};
            Xdpp = obj.parameter{:,'Xd_pp'};
            Xqpp = obj.parameter{:,'Xq_pp'};

            
            syms del e
            Iq = real(I)*cos(del)+imag(I)*sin(del); %Iabs*cos(delta-Iang)
            Id = real(I)*sin(del)-imag(I)*cos(del); %Iabs*sin(delta-Iang)
            Ed = (Xq-Xqp)*Iq;
            psiq = -Ed-(Xqp-Xls)*Iq;
            psid = e-(Xdp-Xls)*Id;
            for_Id = (Xdpp-Xls)*e/(Xdp-Xls) + (Xdp-Xdpp)*psid/(Xdp-Xls);
            for_Iq = (Xqpp-Xls)*Ed/(Xqp-Xls) - (Xqp-Xqpp)*psiq/(Xqp-Xls);
            Vq = -Xdpp*Id + for_Id;
            Vd = Xqpp*Iq+for_Iq;
            
            eq1 = P-Vq*(Iq)-Vd*(Id) == 0;
            eq2 = Q-Vq*(Id)+Vd*(Iq) == 0;
            eq = [eq1;eq2];
            S = solve(eq);
            if(S.e(1)>0);delta = double(S.del(1));Eq = double(S.e(1));
            else; delta = double(S.del(2)); Eq = double(S.e(2));end

            Iq = real(I)*cos(delta)+imag(I)*sin(delta); %Vabs*cos(delta-Vang)
            Id = real(I)*sin(delta)-imag(I)*cos(delta); %Vabs*sin(delta-Vang)

            Ed = (Xq-Xqp)*Iq;
            psiq = -Ed-(Xqp-Xls)*Iq;
            psid = Eq-(Xdp-Xls)*Id;
            
            Vfd = Eq + (Xd-Xdp)*Id;

            
            %{
            delta = Vangle + atan(P/(Q+Vabs^2/Xq));
            Eqnum = P^2*Xdp*Xq + Q^2*Xdp*Xq + Vabs^2*Q*Xq + Vabs^2*Q*Xdp + Vabs^4;
            Eqden = Vabs*sqrt(P^2*Xq^2 + Q^2*Xq^2 + 2*Vabs^2*Q*Xq + Vabs^4);
            Eq = Eqnum/Eqden;
            Ednum = (Xq-Xqp)*Vabs*P;
            Edden = sqrt(P^2*Xq^2 + Q^2*Xq^2 + 2*Vabs^2*Q*Xq +Vabs^4);
            Ed = Ednum/Edden;
            Iq = real(I)*cos(delta)+imag(I)*sin(delta); %Vabs*cos(delta-Vang)
            Id = real(I)*sin(delta)-imag(I)*cos(delta); %Vabs*sin(delta-Vang)

            psiq = Ed-(Xqp-Xls)*Iq;
            psid = Eq-(Xdp-Xls)*Id;
            Vfd = Eq + (Xd-Xdp)*Id;
            %}
            x_avr = obj.avr.initialize(Vfd, Vabs);
            x_gov = obj.governor.initialize(P);
            x_pss = obj.pss.initialize();
            x_st = [delta; 0; Eq; Ed; psiq; psid; x_avr; x_gov; x_pss];
            
            obj.alpha_st = [P; Vfd; Vabs];
            obj.x_equilibrium = x_st;
            
            obj.V_equilibrium = V;
            obj.I_equilibrium = I;
            obj.x_st = x_st;
            obj.V_st = tools.complex2vec(V);
            obj.I_st = tools.complex2vec(I);
            obj.set_linear_matrix(x_st, tools.complex2vec(V));
        end

        function [ret, sys_fb, sys_V] = get_sys(obj, x_st, Vst)
            if nargin < 2 || isempty(x_st)
                x_st = obj.x_st;
            end
            if nargin < 3 || isempty(Vst)
                Vst = obj.V_st;
            end
            omega0 = obj.omega0;
            % パラメータを追加
            Xq = obj.parameter_vec(1);
            Xqp = obj.parameter_vec(2);
            Xqpp = obj.parameter_vec(3);
            Xd = obj.parameter_vec(4);
            Xdp = obj.parameter_vec(5);
            Xdpp = obj.parameter_vec(6);
            Xls = obj.parameter_vec(7);
            Tdo = obj.parameter_vec(8);
            Tqo = obj.parameter_vec(9);
            TTdo = obj.parameter_vec(10);
            TTqo = obj.parameter_vec(11);
            M = obj.parameter_vec(12);
            d = obj.parameter_vec(13);
            
            % 表記を変更：Vabscos -> Vq, Vabssin -> Vd

            % 動揺方程式の状態変数を４つに変更
            % x1 = delta
            % x2 = omega
            % x3 = Eq
            % x4 = Ed
            % x5 = psiq
            % x6 = psid
            A_swing = [0 omega0 0 0 0 0;
                0 -d/M 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0
                0 0 0 0 0 0
                0 0 0 0 0 0];
            % u1 = Pmech
            % u2 = Vfd
            % u3 = Pout
            % u4 = Efd
            % u5 = Efq
            % u6 = psiq_
            % u7 = psid_
            B_swing = [0, 0, 0, 0, 0, 0, 0;
                1/M, 0, -1/M, 0, 0, 0, 0;
                0, 1/Tdo, 0, -1/Tdo, 0, 0, 0;
                0, 0, 0, 0, -1/Tqo, 0, 0;
                0, 0, 0, 0, 0, 1/TTqo, 0;
                0, 0, 0, 0, 0, 0, 1/TTdo];
            % y = [delta, omega, Eq, Ed, psiq, psid]
            C_swing = eye(6);
            sys_swing = ss(A_swing, B_swing, C_swing, 0);
            OutputGroup = struct();
            OutputGroup.delta = 1;
            OutputGroup.omega = 2;
            OutputGroup.Eq = 3;
            OutputGroup.Ed = 4;
            OutputGroup.psiq = 5;
            OutputGroup.psid = 6;
            sys_swing.OutputGroup = OutputGroup;
            InputGroup = struct();
            InputGroup.Pmech = 1;
            InputGroup.Vfd = 2;
            InputGroup.Pout = 3;
            InputGroup.Efd_swing = 4;
            InputGroup.Efq_swing = 5;
            InputGroup.psiq_ = 6;
            InputGroup.psid_ = 7;
            sys_swing.InputGroup = InputGroup;
            
            % テイラー展開を用いるので，ここからは平衡状態
            delta = x_st(1);
            Eq = x_st(3);
            Ed = x_st(4);
            psiq = x_st(5);
            psid = x_st(6);
            
            syms s_delta s_Eq s_Ed s_psiq s_psid s_Vreal s_Vimag
            Vq = s_Vreal*cos(s_delta)+s_Vimag*sin(s_delta);%Vabs*cos(delta-Vang)
            Vd = s_Vreal*sin(s_delta)-s_Vimag*cos(s_delta);%Vabs*sin(delta-Vang)

            for_Id = (Xdpp-Xls)*s_Eq/(Xdp-Xls) + (Xdp-Xdpp)*s_psid/(Xdp-Xls);
            for_Iq = -(Xqpp-Xls)*s_Ed/(Xqp-Xls) + (Xqp-Xqpp)*s_psiq/(Xqp-Xls);
            Id = (for_Id-Vq)/Xdpp;
            Iq = (for_Iq+Vd)/Xqpp;
            
            % |I|cosI, |I|sinIを逆算
            Ir = Id*sin(s_delta)+Iq*cos(s_delta);
            Ii = -Id*cos(s_delta)+Iq*sin(s_delta);
            I = [Ir;Ii]; %1%2

            P = [s_Vreal, s_Vimag]*I; %3
            Efd = s_Eq + (Xd-Xdp)*(Id- ((Xdp-Xdpp)/(Xdp-Xls)^2)*(s_psid+(Xdp-Xls)*Id-s_Eq));%4
            Efq = s_Ed - (Xq-Xqp)*(Iq- ((Xqp-Xqpp)/(Xqp-Xls)^2)*(s_psiq+(Xqp-Xls)*Iq+s_Ed));%5
            psiq_ = -s_psiq-s_Ed-(Xqp-Xls)*Iq; %6
            psid_ = -s_psid+s_Eq-(Xdp-Xls)*Id; %7

            J = jacobian([P, Efd, Efq, psiq_, psid_, I'], ...
                      [s_delta, s_Eq, s_Ed, s_psiq,s_psid, s_Vreal, s_Vimag]);
            sys_fb = ss(eval(subs(J,[s_delta, s_Eq, s_Ed, s_psiq, s_psid, s_Vreal, s_Vimag],[delta, Eq, Ed, psiq, psid, Vst(1), Vst(2)])));
            InputGroup = struct();
            InputGroup.delta = 1;
            InputGroup.Eq = 2;
            InputGroup.Ed = 3;
            InputGroup.psiq = 4;
            InputGroup.psid = 5;
            InputGroup.V = 6:7;
            sys_fb.InputGroup = InputGroup;
            OutputGroup = struct();
            OutputGroup.P = 1;
            OutputGroup.Efd = 2;
            OutputGroup.Efq = 3;
            OutputGroup.psiq_ = 4;
            OutputGroup.psid_ = 5;
            OutputGroup.I = 6:7;
            sys_fb.OutputGroup = OutputGroup; %ok
            
            % これも平衡点
            Vabs = norm(Vst);
            
            sys_V = ss([eye(2); Vst'/Vabs]);
            sys_V.InputGroup.Vin = 1:2;
            OutputGroup = struct();
            OutputGroup.V = 1:2;
            OutputGroup.Vabs = 3;
            sys_V.OutputGroup = OutputGroup; %ok
            
            sys_avr = obj.avr.get_sys(); %avrのシステムを取得
            sys_pss = obj.pss.get_sys(); %pssのシステムを取得
            sys_gov = obj.governor.get_sys(); %governorのシステムを取得
            G = blkdiag(sys_swing, sys_fb, sys_V, sys_avr, -sys_pss, sys_gov);
            ig = G.InputGroup;
            og = G.OutputGroup;
            % Efqを追加、Eq,Edに変更、ig.Vinとog.Vをつなげた
            feedin = [ig.Pout, ig.Efd, ig.Efd_swing, ig.Efq_swing, ig.delta, ig.Eq, ig.Ed, ig.psiq, ig.psid, ig.psiq_, ig.psid_, ig.V, ig.Vabs, ig.Vfd, ig.u_avr, ig.omega, ig.omega_governor, ig.Pmech];
            feedout = [og.P, og.Efd, og.Efd,         og.Efq, og.delta,       og.Eq, og.Ed, og.psiq, og.psid, og.psiq_, og.psid_,og.V, og.Vabs, og.Vfd, og.v_pss, og.omega, og.omega, og.Pmech];
            I = ss(eye(numel(feedin)));
            
            ret = feedback(G, I, feedin, feedout, 1);
        end
    end
end


