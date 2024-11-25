function [Mt,Mp,Me] = my_solve_buxton_numerical(k,TA,tau,R1p,R1e,f,M0,t_span,heaviside_approx)
%k is steepness of the heaviside approximation
%heaviside_approx (boolean) = 0 means just use the exact heaviside function
    H_approx=@(k,t1,t0)(1./(1+exp(-k.*(t1-t0))));
       
    [t,y]=ode113(@(t,y) diff_e(t,y,k,TA,tau,R1p,R1e,f,M0,H_approx,heaviside_approx),t_span,0);
    
    Mt=y;
    Mp=0;
    Me=y;

end

function dydt=diff_e(t,y,k,TA,tau,R1p,R1e,f,M0,H_approx,heaviside_approx)
    if heaviside_approx==1
        ma=H_approx(k,t,TA)*H_approx(k,TA+tau,t).*2.*M0.*exp(-R1p.*TA); %approximately plug-shaped bolus
    else
        ma=heaviside(t-TA).*heaviside(TA+tau - t).*2.*M0.*exp(-R1p.*TA);
    end
        dydt=f*ma-R1e*y; 
end
    