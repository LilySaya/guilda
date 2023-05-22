
Xq = 0.66;
Xqp = 0.05;
Vangle = 0.4;
Vabs = 2.8;
Xdp = 0.02;
P = 2.5;
Q = 3.4;
Iabs = 4;
Iangle = 2.3;
Xdpp = 0.01;
Xqpp = 0.01;
I = 2.3+1.4i;
Xls = 0.01;

            
            syms del e 
            Iq = real(I)*cos(del)+imag(I)*sin(del); %Vabs*cos(delta-Vang)
            Id = real(I)*sin(del)-imag(I)*cos(del); %Vabs*sin(delta-Vang)
            Ed = (Xq-Xqp)*Iq;
            psiq = Ed-(Xqp-Xls)*Iq;
            psid = e-(Xdp-Xls)*Id;

            for_Id = (Xdpp-Xls)*e/(Xdp-Xls) + (Xdp-Xdpp)*psid/(Xdp-Xls);
            for_Iq = -(Xqpp-Xls)*Ed/(Xqp-Xls) + (Xqp-Xqpp)*psiq/(Xqp-Xls);
            Vq = -Xdpp*Id + for_Id
            Vd = Xqpp*Iq-for_Iq

            eq1 = P-Vq*(Iq)-Vd*(Id) == 0;
            eq2 = Q-Vq*(Id)+Vd*(Iq) == 0;
            eq = [eq1;eq2];
            S = solve(eq,[e del])
            S.e
            S.del
            if(S.e(1)>0);delta = double(S.del(1));Eq = double(S.e(1));
            else; delta = double(S.del(2)); Eq = double(S.e(2));end
            %{
            syms del e
            Iq = real(I)*cos(del)+imag(I)*sin(del); %Vabs*cos(delta-Vang)
            Id = real(I)*sin(del)-imag(I)*cos(del); %Vabs*sin(delta-Vang)
            Ed = (Xq - Xqp )*Iq;
            Vd = Ed + Xqp*Iq;
            Vq = e - Xdp*Id;
            eq1 = P-Vq*(Iq)-Vd*(Id) == 0;
            eq2 = Q-Vq*(Id)+Vd*(Iq) == 0;
            eq = [eq1;eq2];
            S = solve(eq);
            if(S.e(1)>0);delta = double(S.del(1));Eq = double(S.e(1));
            else; delta = double(S.del(2)); Eq = double(S.e(2));end
            delta
            Eq
            %}
