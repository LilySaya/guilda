classdef controller_retrofit_LQR_PARK <  controller
    
    properties(Access=private)
       x_avr
       x_pss
       x_gov
    end
    
    properties
        C
        A
        Bw
        Bv
        K
        nx
        Xd
        Xdp
        Xdpp
        Xq
        Xqp
        Xqpp
        Xls
        x0
        V0
        I0
        sys_fb
        sys_design
        avr
        pss
        governor
        Vfd0
        Vabs0
        Efd0
        Efq0
        psiq_0
        psid_0
        Vabscos0
        Eq0
        Ed0
        psiq0
        psid0
        delta0
        x_avr0
    end
    
    methods
        function obj = controller_retrofit_LQR_PARK(net, idx, Q, R, model, model_agc)
            obj@controller(idx, idx);
            if nargin < 5
                model = [];
            end
            if nargin < 6
                model_agc = [];
            end
            obj.avr = net.a_bus{idx}.component.avr;
            obj.pss = net.a_bus{idx}.component.pss;
            obj.governor = net.a_bus{idx}.component.governor;
            n_avr = obj.avr.get_nx;
            n_pss = obj.pss.get_nx;
            n_gov = obj.governor.get_nx;
            nx = 6;
            
            obj.x_avr = @(x) x(nx+(1:n_avr));
            obj.x_pss = @(x) x(nx+n_avr+(1:n_pss));
            obj.x_gov = @(x) x(nx+n_avr+n_pss+(1:n_gov));
            
            if isempty(model)
               model = ss(zeros(2, 5));
               model.InputGroup.Eq_m = 1;
               model.InputGroup.Ed_m = 2;
               model.InputGroup.psiq_m = 3;
               model.InputGroup.psid_m = 4;
               model.InputGroup.delta_m = 5;
               model.OutputGroup.V_m = 1:2;
            end
            
            if isempty(model_agc)
               model_agc = ss(zeros(1, 2));
               model_agc.InputGroup.omega_agc = 1;
               model_agc.InputGroup.delta_agc = 2;
               model_agc.OutputGroup.u_agc = 1;
            end
            
            sys = net.a_bus{idx}.component.get_sys();
            
            sys_ = sys;
            sys_cat = blkdiag(sys, model, model_agc);
            feedout = [sys_cat.OutputGroup.delta, sys_cat.OutputGroup.Eq, sys_cat.OutputGroup.Ed, ...
                sys_cat.OutputGroup.psiq, sys_cat.OutputGroup.psid, sys_cat.OutputGroup.V_m,...
                sys_cat.OutputGroup.u_agc, sys_cat.OutputGroup.omega, sys_cat.OutputGroup.delta];
            feedin = [sys_cat.InputGroup.delta_m, sys_cat.InputGroup.Eq_m, sys_cat.InputGroup.Ed_m,...
                sys_cat.InputGroup.psiq_m, sys_cat.InputGroup.psid_m, sys_cat.InputGroup.Vin,...
                sys_cat.InputGroup.u_governor, sys_cat.InputGroup.omega_agc, sys_cat.InputGroup.delta_agc];
            %sys_cat is G+Gbar(0ony), eye(numel(feedin)) is Gbar
            %feedout is input to Gbar (w)
            %feedin is output from Gbar (v) needs delta
            sys = feedback(sys_cat, eye(numel(feedin)), feedin, feedout, 1);
            
            obj.sys_design = sys;
            
            %A would not be affected whether there is input from Gbar or
            %not
            
            
            [A, B,  C, ~] = ssdata(obj.sys_design('I', {'u_avr'}));
            %[~, N, ~, M] = ssdata(sys('I', {'Vin'}));
            %Pout is P here
            [~, Br, ~, Dr] = ssdata(sys('I', {'Pout', 'Efd_swing', 'Efq_swing','psiq_','psid_','Vabs', 'Vfd' 'u_avr',  'u_governor'}));
            [~, Bw, ~, Dw] = ssdata(sys('I', {'delta', 'delta_m', 'delta_agc', 'omega_agc', 'Eq','Eq_m','Ed','Ed_m','psiq_m','psiq','psid_m','psid', 'Vfd'}));
            
            L = Br;
            obj.C = C;
            obj.A = A;
            obj.Bv = L;
            %sum them, but keep them in columns
            obj.Bw = [sum(Bw(:, [1:3]), 2), Bw(:, 4), sum(Bw(:, [5:6]), 2), sum(Bw(:, [7:8]), 2),sum(Bw(:, [9:10]), 2),sum(Bw(:, [11:12]), 2),Bw(:, 13:end)];

            Q_ = zeros(size(A));
            Q_(1:size(Q, 1), 1:size(Q, 2)) = Q;
            
            if isinf(R)
                obj.K = zeros(1, size(A, 1));
            else
                %No information taken from v, so this K is Khat
                %rectifier should be written somewhr
                obj.K = lqr(A, B(:, 1), Q_, R);
                
            end
            
            obj.nx = size(A, 1);
            obj.Xd = net.a_bus{idx}.component.parameter{:, 'Xd'};
            obj.Xdp = net.a_bus{idx}.component.parameter{:, 'Xd_prime'};
            obj.Xdpp = net.a_bus{idx}.component.parameter{:, 'Xd_pp'};
            obj.Xq = net.a_bus{idx}.component.parameter{:, 'Xq'};
            obj.Xqp = net.a_bus{idx}.component.parameter{:, 'Xq_prime'};
            obj.Xqpp = net.a_bus{idx}.component.parameter{:, 'Xq_pp'};
            obj.Xls = net.a_bus{idx}.component.parameter{:, 'X_ls'};
            obj.x0 = net.a_bus{idx}.component.x_equilibrium;
            obj.V0 = tools.complex2vec(net.a_bus{idx}.component.V_equilibrium);
            obj.I0 = tools.complex2vec(net.a_bus{idx}.component.I_equilibrium);

            %not used in this class
            obj.sys_fb = ss((A-B(:, 1)*obj.K), B(:, 1), [eye(size(A)); -obj.K], 0);

            Vabs0 = norm(obj.V0);
            Xd = obj.Xd;
            Xdp = obj.Xdp;
            Xq = obj.Xq;
            Xqp = obj.Xqp;
            Xdpp = obj.Xdpp;
            Xqpp = obj.Xqpp;
            Xls = obj.Xls;
            delta0 = obj.x0(1); 
            Eq0 = obj.x0(3);
            Ed0 = obj.x0(4);
            psiq0 = obj.x0(5);
            psid0 = obj.x0(6);
            
            Iq0 = obj.I0(1)*cos(delta0)+obj.I0(2)*sin(delta0); %Vabs*cos(delta-Vang)
            Id0 = obj.I0(1)*sin(delta0)-obj.I0(2)*cos(delta0); %Vabs*sin(delta-Vang)

            Efd0 = Eq0 + (Xd-Xdp)*(Id0- ((Xdp-Xdpp)/(Xdp-Xls)^2)*(psid0+(Xdp-Xls)*Id0-Eq0));
            Efq0 = Ed0 - (Xq-Xqp)*(Iq0- ((Xqp-Xqpp)/(Xqp-Xls)^2)*(psiq0+(Xqp-Xls)*Iq0+Ed0));
            [~, obj.Vfd0] = obj.avr.get_Vfd(obj.x_avr(obj.x0), Vabs0, Efd0, 0);
            psiq_0 = -psiq0-Ed0-(Xqp-Xls)*Iq0;
            psid_0 = -psid0+Eq0-(Xdp-Xls)*Id0;

            obj.x_avr0 = obj.x_avr(obj.x0);
            obj.Vabs0 = Vabs0;
            obj.Efq0 = Efq0;
            obj.Efd0 = Efd0;
            obj.psiq_0 = psiq_0;
            obj.psid_0 = psid_0;
            obj.Eq0 = Eq0;
            obj.Ed0 = Ed0;
            obj.psiq0 = psiq0;
            obj.psid0 = psid0;
            obj.delta0 = delta0;
        end
        
        function nx = get_nx(obj)
            nx = obj.nx;
        end
        
        function nu = get_nu(obj)
            nu = 2;
        end
        
        
        function [dx, u] = get_dx_u(obj, t, x, X, V, I, U)
            
            %x: states of component+Gbar 全システムの内部状態
            %X: states of component 機器の内部状態
            %U(1) should be Vpss and U(2) should be Pmech
            u = zeros(2, 1);
            %x0 has 7 attributes and component states
            x1 = x(1:numel(obj.x0));
            x2 = x(numel(obj.x0)+1:end);
            %X{1} has n rows, 7 columns, n depends on time step
            %Why X{1} and not index number of component?
            %Component state - equilibrium - system only affected by environ state
            %C is identity matrix
            u(1) = -obj.K*[( X{1}-obj.x0-x1); -x2];
            
            delta = X{1}(1);
            omega = X{1}(2);
            Eq = X{1}(3);
            Ed = X{1}(4);
            psiq = X{1}(5);
            psid = X{1}(6);
            x_pss = obj.x_pss(X{1});
            x_avr = obj.x_avr(X{1});

            Xd = obj.Xd;
            Xdp = obj.Xdp;
            Xq = obj.Xq;
            Xqp = obj.Xqp;
            Xdpp = obj.Xdpp;
            Xqpp = obj.Xqpp;
            Xls = obj.Xls;
            
            Vabs = norm(V);
            
            P = I'*V - obj.I0'*obj.V0;
            Iq = I(1)*cos(delta)+I(2)*sin(delta) ;%Vabs*cos(delta-Vang)
            Id = I(1)*sin(delta)-I(2)*cos(delta); %Vabs*sin(delta-Vang)

            Efd = Eq + (obj.Xd-obj.Xdp)*(Id- ((obj.Xdp-obj.Xdpp)/(obj.Xdp-obj.Xls)^2)*(psid+(obj.Xdp-obj.Xls)*Id-Eq));
            Efq = Ed - (obj.Xq-obj.Xqp)*(Iq- ((obj.Xqp-obj.Xqpp)/(obj.Xqp-obj.Xls)^2)*(psiq+(obj.Xqp-obj.Xls)*Iq+Ed));
            psiq_ = -psiq-Ed-(obj.Xqp-obj.Xls)*Iq;
            psid_ = -psid+Eq-(obj.Xdp-obj.Xls)*Id;

            [~, vpss] = obj.pss.get_u(x_pss, omega);
            %avr.get_Vfd(V_tr,Vabs,Efd,Vpss)
            %U is input to generator, samll u is output from controller
            [~, Vfd, Vap] = obj.avr.get_Vfd(x_avr, Vabs, Efd, u(1) + U{1}(1)-vpss);
            
            %'Pout', 'Efd_swing', 'Efq_swing','psiq_','psid_','Vabs', 'Vfd' 'u_avr',  'u_governor'
            %'delta', 'delta_m', 'delta_agc', 'omega_agc', 'Eq','Ed','E_m','psiq','psid', 'Vfd'
            %IV, Efd, Vabs, Vfd 'Pout', 'Efd_swing', 'Vabs', 'Vfd' 'u_avr',  'u_governor'
            %delta, omega, E, Vap 'delta', 'delta_m', 'delta_agc', 'omega_agc', 'E', 'E_m', 'Vfd'
            %System only affected by environment like Gyv, Gyw
            %Bv 9inputs-governor
            %Pmech = obj.governor.get_P(obj, x_gov, omega, U{1}(2));

            dx = obj.A*x + obj.Bv*[P; Efd-obj.Efd0; Efq-obj.Efq0; psiq_-obj.psiq_0; psid_-obj.psid_0; Vabs-obj.Vabs0; Vfd-obj.Vfd0; U{1}] ...
                - obj.Bw*[delta-obj.delta0; omega; Eq-obj.Eq0; Ed-obj.Ed0; psiq-obj.psiq0; psid-obj.psid0; Vap-obj.Vfd0];
        end
    end
end
