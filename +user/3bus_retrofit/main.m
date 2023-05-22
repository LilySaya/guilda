%simulation conditions
%time = [0,100];
%option.fault = {{[10 10.01], 1}};
%parameters - machinery.xlsx 
%machinery_same.xlsx for very small TTqo, TTdo

%No control
%no_control

%AGC only
%AGC_only

%AGC and PSS and AVR
%AGC_PSS_AVR

%AGC and PSS and AVR and retrofit

%{
axis1_lqr
hold on;
axis2_lqr
hold on;
park_lqr
legend('1axis','2axis','park');
%}


%1axis
%axis1_comparison

%2axis
%axis2_comparison

%park
%park_comparison

plot(out.t,out.X{1}(:,2));
hold on;
plot(out2.t,out2.X{1}(:,2));
hold on;
plot(out3.t,out3.X{1}(:,2));
legend('1axis','2axis','park');


