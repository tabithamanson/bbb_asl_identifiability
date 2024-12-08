% Fitting simulations. Series 2CXM, multi-echo scan.

n=500; % no. of instances
default_min=eps;
k=100; % Heaviside steepness
heaviside_approx=0;
tau=1.000; % BL
reps=2;
TA_and_f_uncertainty=0;
use_identifiable_version=0; %THIS HAS BUGS instead change the lower and upper bounds on the objective function

R1p_uncertainty=1; %
R2e_uncertainty=1; % need to change these in the actual bounds, too

fit_R2e=0; % NOT IMPLEMENTED
A=500; %exponential constant term for monoexponential R2e fit (not implemented)

use_series_model=1;

noise=1;

tic
%Multi-echo t_array
TI=[1.100 2.100 3.100];
tes=[20.8, 62.5, 104.2, 145.8, 187.5, 229.2, 270.9]./1000;
t_array=ones(1,21);
t_array=[TI(1).*ones(1,7) TI(2).*ones(1,7) TI(3).*ones(1,7)]+[tes tes tes];
t_span=t_array;

run literature_vals.m

%pre-allocations:
f_GT=ones(1,n).*nan;
TA_GT=ones(1,n).*nan;
R1p_GT=ones(1,n).*nan;
R1e_GT=ones(1,n).*nan;
R2p_GT=ones(1,n).*nan;
R2e_GT=ones(1,n).*nan;
kw_GT=ones(1,n).*nan;
d_GT=ones(1,n).*nan;
x_GT=ones(n,5).*nan;
data=ones(n,42).*nan;
data_R2=ones(n,7).*nan;
TA_est=ones(1,n).*nan;
f_est=ones(1,n).*nan;
R1p_est=ones(1,n).*nan;
R2e_est=ones(1,n).*nan;

x0=[TA+1/kw R1p R1e R2p R2e];

x_fit=nan*ones(n,length(x0));
RESNORM=nan*ones(n,1);
RESIDUAL=nan*ones(n,1);
EXITFLAG=nan*ones(n,1);
for i=1:n
% Generate GT parameter values (5 parameters: f, R1p, R1e, kw, TA) in literature range (normal dist.). 
    T1p=1/R1p;
    T1e=1/R1e;
    T2p=1/R2p;
    T2e=1/R2e;
    f_GT(i) = max(default_min,normrnd(f,f*0.2));
    TA_GT(i) = max(default_min,normrnd(TA,TA*0.15));
    R1p_GT(i)=1/max(default_min,normrnd(T1p,T1p.*0.05));
    R1e_GT(i)=1/max(default_min,normrnd(T1e,T1e.*0.05));
    R2p_GT(i)=1/max(default_min,normrnd(T2p,T2p.*0.1));
    R2e_GT(i)=1/max(default_min,normrnd(T2e,T2e.*0.2));
    kw_GT(i)=max(default_min,rand*500/60);
    Texch_GT(i)=1./kw_GT(i);%max(default_min,rand*1.10);
    d_GT(i)=Texch_GT(i) + TA_GT(i);

    % store x_GT for this instance
    x_GT(i,:)=[d_GT(i) R1p_GT(i) R1e_GT(i) R2p_GT(i) R2e_GT(i)];

    % Simulated clean signal from this "voxel"
 
            data(i,:) =repmat([my_solve_series_2CXM_T2_numerical(k,TA_GT(i),tau,R1p_GT(i),R1e_GT(i),R2p_GT(i),R2e_GT(i),f_GT(i),d_GT(i),M0,t_span(1:7),1.100,0),...
                my_solve_series_2CXM_T2_numerical(k,TA_GT(i),tau,R1p_GT(i),R1e_GT(i),R2p_GT(i),R2e_GT(i),f_GT(i),d_GT(i),M0,t_span(8:14),2.100,0),...
                my_solve_series_2CXM_T2_numerical(k,TA_GT(i),tau,R1p_GT(i),R1e_GT(i),R2p_GT(i),R2e_GT(i),f_GT(i),d_GT(i),M0,t_span(15:21),3.100,0)],1,reps);


if noise==1
    for k=1:size(data,2)
        %data(i,k)=max(eps,normrnd(data(i,k),sqrt(max(eps,data(i,k).*0.95)))); %if voxelwise
        data(i,k)=max(eps,normrnd(data(i,k),sqrt(max(eps,data(i,k).^1.3.*exp(-3.4))))); % if whole_GM 
    end
end

% Fix R1p & R1e, TA at GT 
if TA_and_f_uncertainty
    TA_est(i)=max(default_min,normrnd(TA_GT(i),0.14*1.000));
    f_est(i)=max(default_min,normrnd(f_GT(i),9.5/(60*100*1.000))); %12, or 9.5 if using single long PLD to estimate
else
    TA_est(i)=TA_GT(i);
    f_est(i)=f_GT(i);
end

if R1p_uncertainty
    R1p_est(i)=R1p;
else
    R1p_est(i)=R1p_GT(i);
end

if R2e_uncertainty
    R2e_est(i)=R2e;
else
    R2e_est(i)=R2e_GT(i);
end

% Fit R2e from L-weighted signal
if fit_R2e==1
objective_function_R2e=@(z) z(2).*exp(-z(1).*tes) + z(3) - data_R2e(i,:);
z0=[R2e 1 0];
opts = optimset('Display','off','Algorithm','trust-region-reflective','FinDiffType','central');
z_fit(i,:)=lsqnonlin(objective_function_R2e,z0,[],[],opts);
else
    z_fit(i,:)=R2e_GT(i);
end
%weights=ones(1,length(data(i,:))).*(1./(0.95.*max(data(i,:),eps)));
% Fit kw and f to generated data (store in n x 5 matrix pfit_vals)

    
     %objective_function = @(x) repmat([my_solve_series_2CXM_T2_numerical(k,TA_est(i),tau,x(2),x(3),x(4),x(5),f_est(i),x(1),M0,t_span(1:7),1.1,0),...
      %              my_solve_series_2CXM_T2_numerical(k,TA_est(i),tau,x(2),x(3),x(4),x(5),f_est(i),x(1),M0,t_span(8:14),2.100,0),...
       %             my_solve_series_2CXM_T2_numerical(k,TA_est(i),tau,x(2),x(3),x(4),x(5),f_est(i),x(1),M0,t_span(15:21),3.100,0)],1,reps)...
        %            -data(i,:);
    objective_function_log=@(x) log(repmat([my_solve_series_2CXM_T2_numerical(k,TA_est(i),tau,x(2),x(3),x(4),x(5),f_est(i),x(1),M0,t_span(1:7),1.100,0),...
                    my_solve_series_2CXM_T2_numerical(k,TA_est(i),tau,x(2),x(3),x(4),x(5),f_est(i),x(1),M0,t_span(8:14),2.100,0),...
                    my_solve_series_2CXM_T2_numerical(k,TA_est(i),tau,x(2),x(3),x(4),x(5),f_est(i),x(1),M0,t_span(15:21),3.100,0)],1,reps) +1)...
                    -log(data(i,:)+1);



    lb=[TA_est(i) R1p R1e R2p R2e_GT(i)];
    ub=[+Inf R1p R1e R2p R2e_GT(i)];
    
    opts = optimset('Display','off','Algorithm','trust-region-reflective','FinDiffType','central');%,'TolFun',1e-9);

 [x_fit(i,:),RESNORM(i),~,EXITFLAG(i)]=lsqnonlin(objective_function_log,x0,lb,ub,opts);


if rem(i,10)==0
disp(['i = ' num2str(i)])
end
end
toc

% figure
% plot(x_GT(:,1)*60*100*1000,x_fit(:,1)*60*100*1000,'o')
% xlabel('f ground-truth')
% ylabel('fitted f')
if use_identifiable_version
    ARE_d=100.*median(abs((x_GT(:,1))-(x_fit(:,1)))./(x_GT(:,1)))
    ARE_R1e=100.*median(abs(x_GT(:,2)-x_fit(:,2))./x_GT(:,2))
    ARE_R2p=100.*median(abs(x_GT(:,3)-x_fit(:,3))./x_GT(:,3))
    %ARE_R2e=100.*median(abs(x_GT(:,4)-x_fit(:,4))./x_GT(:,4))%(R2e_GT - z_fit(:,1).')./R2e_GT)

    figure
    plot(1./x_GT(:,1),1./x_fit(:,1),'x')
    xlabel('1/T_{exch} ground-truth (s)')
    ylabel('fitted 1/T_{exch} (s)')
    
    figure
    plot(1./x_GT(:,2),1./x_fit(:,2),'x')
    xlabel('T1_e ground-truth [s]')
    ylabel('fitted T1_e [s]')
    
    figure
    plot(1./x_GT(:,3),1./x_fit(:,3),'x')
    xlabel('T2_p ground-truth [s]')
    ylabel('fitted T2_p [s]')
    
else

    figure('Units','centimeters','Position',[0, 0, 8.78, 6.6])
    plot(x_GT(:,1),x_fit(:,1),'.')
    xlabel('\delta_t ground-truth (s)')
    ylabel('fitted \delta_t (s)')
    ylim([0 1000])
   

    Texch_fit=(x_fit(:,1)-TA_est.').';

    figure('Units','centimeters','Position',[0, 0, 8.78, 6.6])
    plot(60./Texch_GT,60./Texch_fit,'.')
    hold on
    plot([0:1:500],[0:1:500],'k-')
    xlabel('1/T_{exch} ground-truth (min^{1})')
    ylabel('fitted 1/T_{exch} (min^{1})')
    ylim([0 1000])
    xlim([0 500])
    
    kw_over_1000=sum(60./Texch_fit>1000)
    
    figure
    plot(1./x_GT(:,2),1./x_fit(:,2),'x')
    xlabel('T1_b ground-truth (s)')
    ylabel('fitted T1_b (s)')
    
    figure('Units','centimeters','Position',[0, 0, 8.78, 6.6])
    plot(1./x_GT(:,3),1./x_fit(:,3),'x')
    xlabel('T1_t ground-truth (s)')
    ylabel('fitted T1_t (s)')
    
    figure('Units','centimeters','Position',[0, 0, 8.78, 6.6])
    plot(1./x_GT(:,4),1./x_fit(:,4),'x')
    xlabel('T2_b ground-truth (s)')
    ylabel('fitted T2_b (s)')
    
    figure('Units','centimeters','Position',[0, 0, 8.78, 6.6])
    plot(1./x_GT(:,5),1./x_fit(:,5),'x')
    xlabel('T2_t ground-truth (s)')
    ylabel('fitted T2_t (s)')
    
    %Calculate average relative estimation error  (ARE) for kw and f over all iterations
    % RE_f=(x_GT(:,1)-x_fit(:,1))./x_GT(:,1);
    % ARE_f=100.*(1/n).*(x_GT(:,1)-x_fit(:,1))./x_GT(:,1)
    
    RE_Texch=100.*median(((x_GT(:,1))-(x_fit(:,1)))./(x_GT(:,1)))
    RE_R1p=100.*median((x_GT(:,2)-x_fit(:,2))./x_GT(:,2))
    RE_R1e=100.*median((x_GT(:,3)-x_fit(:,3))./x_GT(:,3))
    RE_R2p=100.*median((x_GT(:,4)-x_fit(:,4))./x_GT(:,4))
    RE_R2e=100.*median((x_GT(:,5)-x_fit(:,5))./x_GT(:,5))
    
    ARE_d=100.*median(abs((x_GT(:,1))-(x_fit(:,1)))./(x_GT(:,1)))
    ARE_Texch=100.*median(abs(60./Texch_GT-60./Texch_fit)./(60./Texch_GT))
    ARE_R1p=100.*median(abs(x_GT(:,2)-x_fit(:,2))./x_GT(:,2))
    ARE_R1e=100.*median(abs(x_GT(:,3)-x_fit(:,3))./x_GT(:,3))
    ARE_R2p=100.*median(abs(x_GT(:,4)-x_fit(:,4))./x_GT(:,4))
    ARE_R2e=100.*median(abs(x_GT(:,5)-x_fit(:,5))./x_GT(:,5))%(R2e_GT - z_fit(:,1).')./R2e_GT)

end

%plot error propagation graph
% figure
% semilogy(100.*(abs(TA_est-TA_GT)./TA_GT),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.');
% xlabel('Absolute Relative Error in T_A %')
% ylabel('Absolute Relative Error in k_w %')
% 
% figure
% semilogy(100.*(abs(f_est-f_GT)./f_GT),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.');
% xlabel('Absolute Relative Error in f %')
% ylabel('Absolute Relative Error in k_w %')
% 
% figure
% semilogy(100.*(abs(R1p_est-R1p_GT)./R1p_GT),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.')
% xlabel('Absolute Relative Error in R1_b %')
% ylabel('Absolute Relative Error in k_w %')

% x_lit=[kw R1p R1e R2p R2e];
% data_lit=objective_function_log(x_lit)
% objective_function(x_lit)m 

figure
plot(100.*(abs(x_fit(:,3)-x_GT(:,3))./x_GT(:,3)),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.')
xlabel('Absolute Relative Error in R1_e %')
ylabel('Absolute Relative Error in d_t %')

figure
plot(100.*(abs(x_fit(:,4)-x_GT(:,3))./x_GT(:,4)),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.')
xlabel('Absolute Relative Error in R2_p %')
ylabel('Absolute Relative Error in d_t %')

% figure
% semilogy(100.*(abs(R2p_est-R2p_GT)./R2p_GT),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.')
% xlabel('Absolute Relative Error in R2_p %')
% ylabel('Absolute Relative Error in k_w %')
% 
% figure
% semilogy(100.*(abs(R2e_est-R2e_GT)./R2e_GT),100.*abs(x_GT(:,1)-x_fit(:,1))./x_GT(:,1),'.')
% xlabel('Absolute Relative Error in R2_e %')
% ylabel('Absolute Relative Error in k_w %')

% figure
% plot(100.*(x0(:,1)-x_GT(:,1))./x_GT(:,1),RE_kw*100,'x')
% xlabel('kw initial point deviation from ground truth (%)')
% ylabel('kw % relative error')

% figure
% plot(100.*(x0(:,2)-x_GT(:,2))./x_GT(:,2),RE_kw*100,'x')
% xlabel('kw initial point deviation from ground truth (%)')
% ylabel('kw % relative error')

%workspace_name=['./series2CXM_sim_fit_noise_' num2str(noise) '_GM_RA_dt'];
%save(workspace_name);
