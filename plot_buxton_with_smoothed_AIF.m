%Plot buxton 1 CM with analytical solution vs numerical heaviside vs
%numerical with heaviside approximation 

tau=0.400; %s
t_span=(0:10:5000)./1000;
t=t_span;
run("literature_vals.m")
k=100;
model_fun=load('my_buxton_analytical.mat');
Mt_analytical=double(subs(model_fun.Mt));
Mt_numerical = my_solve_buxton_numerical(0,TA,tau,R1p,R1e,f,M0,t_span,0);
Mt_numerical_k01 = my_solve_buxton_numerical(100,TA,tau,R1p,R1e,f,M0,t_span,1);
Mt_numerical_k005 = my_solve_buxton_numerical(50,TA,tau,R1p,R1e,f,M0,t_span,1);
Mt_numerical_k0025= my_solve_buxton_numerical(25,TA,tau,R1p,R1e,f,M0,t_span,1);

figure('units','centimeters','position',[0,0,13.02,2*13.02/3])
plot(t,Mt_analytical,'k-',...
    t,Mt_numerical,'c-',...
    t,Mt_numerical_k01,'r-',...
    t,Mt_numerical_k005,'g-',...
    t,Mt_numerical_k0025,'b-','LineWidth', 1)
legend('analytical solution','numerical solution','c=100','c=50','c=25')
xlabel('time (s)')
ylabel("magnetisation (AU)")
xlim([1.300 3])

alw = 1;    % AxesLineWidth
fsz = 12;      % Fontsize
set(gca, 'FontSize', fsz, 'LineWidth', alw)
set(gcf,'Color','w')
