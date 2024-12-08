function [Mt,Mp,Me,Mt0] = my_solve_series_2CXM_T2_numerical(k,TA,tau,R1p,R1e,R2p,R2e,f,d,M0,t_array,TI_array,heaviside_approx)
% Solve series_2CXM numerically (multi-echo scan context)
%k is steepness of the heaviside approximation
%heaviside_approx (boolean) = 0 means just use the exact heaviside function
%NOTE I typically call this with just one TI so multiple TIs in "TI_array" might be buggy 

%TI_array should be given to the nearest 10ms (otherwise adjust TItime
%below)

vp=0.05;

H_approx=@(k,t1,t0)(1./(1+exp(-k.*(t1-t0))));
num_echoes=length(t_array)/length(TI_array);
    for i=1:length(TI_array)
        % Find index of t=TI
        TI0=TI_array(i);
        T1t=(0:10:5100)./1000; % 5100 arbitrarily chosen maximum time (increase if needed). Change increment from 10 to 1 if TI is given to the nearest ms (currently to the nearest 10 ms)
        ind0=find(T1t==TI0);
        % Find longitudinal magn. at t=TI (init. conditions)
        [Mt0,Mp0,Me0]=my_solve_series_2CXM_numerical(k,TA,tau,R1p,R1e,f,d,M0,(0:10:5000)./1000, 1);
        Mt0=Mp0(ind0)+Me0(ind0);
    
        % Solve
        [t,y]=ode113(@(t,y) diff_m(t,y,k,TA,tau,R1p,R1e,R2p,R2e,f,d,vp,M0,H_approx,heaviside_approx),t_array((i-1)*num_echoes+1:(i)*num_echoes),[Mp0(ind0)./vp,Me0(ind0)./(1-vp)]);
    
        Mt((i-1)*num_echoes+1:(i)*num_echoes)=vp.*y(:,1) + (1-vp).*y(:,2);
        Mp((i-1)*num_echoes+1:(i)*num_echoes)=vp.*y(:,1);
        Me((i-1)*num_echoes+1:(i)*num_echoes)=(1-vp).*y(:,2);
    end
end

function dydt=diff_m(t,y,k,TA,tau,R1p,R1e,R2p,R2e,f,d,vp,M0,H_approx,heaviside_approx)
dydt=zeros(2,1);
if heaviside_approx==1
    ma=H_approx(k,t,TA)*H_approx(k,TA+tau,t).*2.*M0.*exp(-R1p.*TA); %approximately plug-shaped bolus
else
    ma=heaviside(t-TA).*heaviside(TA+tau - t).*2.*M0.*exp(-R1p.*TA);
end
ma=0;
 dydt(1)=-(vp.*R2p.*y(1))./vp;
 dydt(2)=-(1-vp).*R2e.*y(2)./(1-vp);
end
    