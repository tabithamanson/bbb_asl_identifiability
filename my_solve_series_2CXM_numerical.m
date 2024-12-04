function [Mt,Mp,Me] = my_solve_series_2CXM_numerical(k,TA,tau,R1p,R1e,f,d,M0,t_span,heaviside_approx)
% Solve series_2CXM numerically (single-echo scan context)
%k is steepness of the heaviside approximation
%heaviside_approx (boolean) = 0 means just use the exact heaviside function
vp=0.05;
H_approx=@(k,t1,t0)(1./(1+exp(-k.*(t1-t0))));

[t,y]=ode113(@(t,y) diff_m(t,y,k,TA,tau,R1p,R1e,f,d,vp,M0,H_approx,heaviside_approx),t_span,[0 0]);

Mt=vp.*y(:,1) + (1-vp).*y(:,2); %total signal
Mp=vp.*y(:,1); %blood signal
Me=(1-vp).*y(:,2); %tissue signal
%t_span=t;
%save('my_buxton_smooth_numerical.mat','Mt','Mp','Me')
end

function dydt=diff_m(t,y,k,TA,tau,R1p,R1e,f,d,vp,M0,H_approx,heaviside_approx)
dydt=zeros(2,1);
%d=Texch + TA;
if heaviside_approx==1
    ma=H_approx(k,t,TA)*H_approx(k,TA+tau,t).*2.*M0.*exp(-R1p.*TA); %approximately plug-shaped bolus
    me_in=H_approx(k,t,d).*2.*M0.*exp(-R1p.*d).*H_approx(k,d+tau,t);
else
    ma=heaviside(t-TA).*heaviside(TA+tau - t).*2.*M0.*exp(-R1p.*TA);
    me_in=heaviside(t-d).*2.*M0.*exp(-R1p.*d).*heaviside(d+tau-t);
end
dydt(1)=((f*ma -R1p.*y(1).*vp)*H_approx(k,d,t) - (f*me_in+y(1)*R1p*vp)*H_approx(k,t,TA+tau))./vp;
dydt(2)=-(1-vp)*R1e*y(2)/(1-vp)+ f*me_in/(1-vp);
end
    