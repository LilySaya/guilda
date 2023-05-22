%{
axis1_hinf
hold on;
axis2_hinf
hold on;
park_hinf
%}

%{
axis1_lqr
axis2_lqr
park_lqr
%}

%{
init_axis1
hold on;
init_axis2
hold on;
init_park
%}

plot(out.t,out.X{1}(:,2));
hold on;
plot(out2.t,out2.X{1}(:,2));
hold on;
plot(out3.t,out3.X{1}(:,2));
legend('1axis','2axis','park');

