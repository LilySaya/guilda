classdef generator_2axis < component % 状態・パラメーターはqを先においている
    %edit
    %forgot but edited +,- signs somewhere
    properties(Access = private)
        parameter_vec
        system_matrix
        x_st
        V_st
        I_st
    end
    
    properties(SetAccess = private)
%         parameter
%         x_equilibrium
%         V_equilibrium
%         I_equilibrium
        avr
        pss
        governor
        alpha_st
        omega0
    end
    
    methods
        function obj = generator_2axis(omega, parameter)
            obj.omega0 = omega;
            if isstruct(parameter)
                parameter = struct2table(parameter);
            end
            % 2軸用のパラメータ名に変更
            obj.parameter = parameter(:, {'Xq', 'Xq_prime', 'Xq_pp','Xd', 'Xd_prime', 'Xd_pp','X_ls','Tdo', 'Tqo', 'TTdo','TTqo','M', 'D'});   % ソートしてるだけ
            obj.parameter_vec = obj.parameter.Variables;
            obj.avr = avr();
            obj.governor = governor();
            obj.pss = pss();
            obj.system_matrix = struct();
        end
        
        function name_tag = get_x_name(obj)
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

            % Id, Iqを定義
            Iq = -(Ed-Vd)/Xqp;
            Id = (Eq-Vq)/Xdp;
            
            % |I|cosI, |I|sinIを逆算
            Ir = Id*sin(delta)+Iq*cos(delta);
            Ii = -Id*cos(delta)+Iq*sin(delta);
            
            con = I - [Ir; Ii];
            
            % Efdの修正とEfqの追加
            Efd = Eq + (Xd-Xdp)*(Id- ((Xdp-Xdpp)/(Xdp-Xls)^2)*(psid+(Xdp-Xls)*Id-Eq));
            Efq = Ed - (Xq-Xqp)*(Iq- ((Xqp-Xqpp)/(Xqp-Xls)^2)*(psiq+(Xqp-Xls)*Iq+Ed));
            
            [dx_pss, v] = obj.pss.get_u(x_pss, omega);
            [dx_avr, Vfd] = obj.avr.get_Vfd(x_avr, Vabs, Efd, u(1)-v);
            [dx_gov, P] = obj.governor.get_P(x_gov, omega, u(2));
            
            % dEをdEqに，dEdの追加
            dEq = (-Efd + Vfd)/Tdo;
            dEd = (-Efq)/Tqo;
            ddelta = omega0 * omega; 
            % PはPmechを指す
            domega = (P - d*omega - Vq*Iq - Vd*Id)/M;

            %dpsid, dpsiqを追加
            dpsiq = (-psiq-Ed-(Xqp-Xls)*Iq)/TTqo;
            dpsid = (-psid+Eq-(Xdp-Xls)*Id)/TTdo;
            % ここで，dEqとdEdの順序を逆にする必要があるかも
            dx = [ddelta; domega; dEq; dEd; dpsiq; dpsid; dx_avr; dx_pss; dx_gov];
        end %ok
        
       
        
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
            
            delta = Vangle + atan(P/(Q+Vabs^2/Xq));
            Eqnum = P^2*Xdp*Xq + Q^2*Xdp*Xq + Vabs^2*Q*Xq + Vabs^2*Q*Xdp + Vabs^4;
            Eqden = Vabs*sqrt(P^2*Xq^2 + Q^2*Xq^2 + 2*Vabs^2*Q*Xq + Vabs^4);
            Eq = Eqnum/Eqden;
            Ednum = (Xq-Xqp)*Vabs*P;
            Edden = sqrt(P^2*Xq^2 + Q^2*Xq^2 + 2*Vabs^2*Q*Xq +Vabs^4);
            Ed = Ednum/Edden;
            

            psiq = (-Ed-(Xqp-Xls)*Iabs*cos(delta-Iangle));
            psid = (Eq-(Xdp-Xls)*Iabs*sin(delta-Iangle));


            Vfd = Eq + (Xd-Xdp)*Iabs*sin(delta-Iangle);

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
            %obj.set_linear_matrix(x_st, tools.complex2vec(V));
        end
    end
end


