function [Mt,Mp,Me] = my_solve_parallel_2CXM_numerical(k,TA,tau,R1p,R1e,f,kw,M0,t_span,heaviside_approx)
% Solve parallel_2CXM numerically (single-echo scan context)
%k is steepness of the heaviside approximation
%heaviside_approx (boolean) = 0 means just use the exact heaviside function
vp=0.05;
H_approx=@(k,t1,t0)(1./(1+exp(-k.*(t1-t0))));
   
[t,y]=ode113(@(t,y) diff_m(t,y,k,TA,tau,R1p,R1e,f,kw,vp,M0,H_approx,heaviside_approx),t_span,[0 0]);
    Mt=vp.*y(:,1) + (1-vp).*y(:,2); %total signal
    Mp=vp.*y(:,1); %blood signal
    Me=(1-vp).*y(:,2); %EES signal
    t=t_span;
end

function dydt=diff_m(t,y,k,TA,tau,R1p,R1e,f,kw,vp,M0,H_approx,heaviside_approx)
dydt=zeros(2,1);
if heaviside_approx==1
    ma=H_approx(k,t,TA)*H_approx(k,TA+tau,t).*2.*M0.*exp(-R1p.*TA); %approximately plug-shaped bolus
else
    ma=heaviside(t-TA).*heaviside(TA+tau - t).*2.*M0.*exp(-R1p.*TA);
end
 dydt(1)=f*ma/vp-R1p*y(1)-kw*y(1);
 dydt(2)=-R1e*y(2)+ vp*kw*y(1)/(1-vp);
end
    