% Structural identifiability analysis using the sensitivity matrix.
% Produces identifiability signatures and sensitivity plots, and lists
% locally non-identifiable parameters.
% Check for convergence first by changing 'm' and assessing "max_change"
% output.

% CHECK before running: 
% m : this controls the minimum step size. Start at e.g. 10 and reduce until the output "max_change" < 0.01 between iterations.
% TI : (t_array and subbed values)
% t_array : (acquisition times ; single-te TI times or multi-te TI+TE times)
% p : free parameters (exclude fixed parameters)
% tau : bolus duration (1000 ms for multi-te and 400 ms for single-te)
% Fixed parameter values in function at end of file
% Uncomment the correct model (2x) see lines marked with %++++++++%
% clear workspace before running 

clear workspace
sympref('TypesetOutput',true);
syms x_dot(t) g(t) y(t) x(t) Y(s)
syms x1p x1e positive real rational
syms M0 f F TA R1p kw tau R1e vp t R2p R2e k positive real rational
syms s positive

alw = 2;    % AxesLineWidth
fsz = 12;      % Fontsize

% make a struct of free parameters
p=[f TA R1p R1e]; % kw];% M0 tau];
m_max=4; 

% initial condition
x0=[0; 0];

vars=[p];

num_p=length(p);
num_x=size(x0,1);

% scale factor big "C"
C=[1 1]; % no scaling

y(t,p)=C*x;

% single-TE scan t_array:
t_array=[500:400:2900]./1000; %s

% Multi-echo scan t_array:
% TI=3100;
% tes=[20.8, 62.5, 104.2, 145.8, 187.5, 229.2, 270.9];
% t_array=ones(1,7);
% t_array(1:7)=TI*ones(1,7)+tes;

num_t=length(t_array);
tau=0.400;
run('literature_vals.m')
c=100; % heaviside approximation steepness
heaviside_approx=1;

% create matrix of central differences about literature values:
syms y_subs(t)
y_subs(t)=subs(y);
p_subs=subs(p);

%left parameter values = 0.99*param vals
%right param vals = 1.01*param vals
change(1)=+Inf; % first change in S entries (from "no value" to something)
change_max(1)=+Inf;
m=0;
p_del=+Inf;
while (sum(p_del<10*eps) < 1 && m<m_max)
    m=m+1;
    del_L=-(1e-2)*1*10^(-(m-1)); % little step to the left in parameter space
    del_R=-del_L; % little step to the right in parameter space
    for k=1:(num_p)
        switch p(k)
            case 'f'
                p_del(k)=(del_R*f);
                f=(1+del_L).*f;
            case 'TA'
                p_del(k)=(del_R*TA);
                TA=(1+del_L).*TA;
            case 'R1p'
                p_del(k)=(del_R*R1p);
                R1p=(1+del_L).*R1p;
            case 'R1e'
                p_del(k)=(del_R*R1e);
                R1e=(1+del_L).*R1e;
            case 'R2p'
                p_del(k)=(del_R*R2p);
                R2p=(1+del_L).*R2p;
            case 'R2e'
                p_del(k)=(del_R*R2e);
                R2e=(1+del_L).*R2e;
            case 'kw'
                p_del(k)=(del_R*kw);
                kw=(1+del_L).*kw;
            otherwise
                disp('error in the switch statement :( ')
        end
        [f,TA,R1p,R1e,R2p,R2e]=fixparameters(f,TA,R1p,R1e,R2p,R2e);
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        [Mt_L,Mp_L,Me_L] = my_solve_buxton_numerical(c,TA,tau,R1p,R1e,f,M0,t_array,heaviside_approx);
        %[Mt_L,Mp_L,Me_L] = my_solve_parallel_2CXM_numerical(c,TA,tau,R1p,R1e,f,kw,M0,t_array,heaviside_approx);
        %[Mt_L,Mp_L,Me_L] = my_solve_series_2CXM_numerical(c,TA,tau,R1p,R1e,f,(1/kw + TA),M0,t_array,heaviside_approx);
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        
        t=t_array;
        M(:,1,k)=Mt_L;
        run('literature_vals.m')
    end

    for k=1:(num_p)
        switch p(k)
            case 'f'
                f=(1+del_R).*f;
            case 'TA'
                TA=(1+del_R).*TA;
            case 'R1p'
                R1p=(1+del_R).*R1p;
            case 'R1e'
                R1e=(1+del_R).*R1e;
            case 'R2p'
                R2p=(1+del_R).*R2p;
            case 'R2e'
                R2e=(1+del_R).*R2e;
            case 'kw'
                kw=(1+del_R).*kw;
            otherwise
                disp('error in the switch statement :( ')
        end
        [f,TA,R1p,R1e,R2p,R2e]=fixparameters(f,TA,R1p,R1e,R2p,R2e);
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        [Mt_R,Mp_R,Me_R] = my_solve_buxton_numerical(c,TA,tau,R1p,R1e,f,M0,t_array,heaviside_approx);
        %[Mt_R,Mp_R,Me_R] = my_solve_parallel_2CXM_numerical(c,TA,tau,R1p,R1e,f,kw,M0,t_array,heaviside_approx);
        %[Mt_R,Mp_R,Me_R] = my_solve_series_2CXM_numerical(c,TA,tau,R1p,R1e,f,(1/kw + TA),M0,t_array,heaviside_approx);
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

        M(:,2,k)=Mt_R;
        t=t_array;
        run('literature_vals.m')

    end

    %if the estimated signal is close to machine zero, round it down to
    %zero:
    M(abs(M)<10*eps)=0;

    for k=1:(num_p)
        Dfull(:,k) = (M(:,2,k) - M(:,1,k))./(2*p_del(k));

        if(p_del(k)<10*eps)
            disp('warning: step size < eps')
        end
    end

    Dfull(abs(Dfull)<10*eps)=0;
    D=Dfull;
    S_sub(:,:,m)=D;

    if m>1
        S_change(:,:,m-1) = abs(S_sub(:,:,m)-S_sub(:,:,m-1))./max(eps,abs(S_sub(:,:,m-1)));
        change(m)=mean(S_change(:,:,m-1),'all');
        change_max(m)=max(S_change(:,:,m-1),[],'all');
    end

    xd(m)=del_R*100;
    [U SIG VT]=svd(S_sub(:,:,m));
    VT(abs(VT)<0.001)=0;

    singular_vals=SIG(1:num_p,1:num_p)*ones(size(SIG(1:num_p,1:num_p),1),1);
    VT_small=[];
    V_small=[];
    for i=2:length(singular_vals)
        if log10(singular_vals(i))<log10(singular_vals(i-1))-3
            VT_small=VT(:,i:end);
            V_small=VT_small.';
            break
        end
    end

    unid_params=[];

    for i=1:size(V_small,1)
        nonzero_V_small=abs(V_small(i,:))>eps;
        unid_params=[unid_params,string(p(nonzero_V_small))];
        unid_params=unique(unid_params);
    end

    disp(["     m    del        max_change      rank    rank_eff    unid_params"])
    disp([string(m)   string(xd(m))      string(change_max(m))    string(rank(S_sub(:,:,m)))   string(size(VT,1)-size(V_small,1))           unid_params])

    SIG(abs(SIG)<=10*eps)=0;
end

%get singular vals from 1:num_p
singular_vals=SIG(1:num_p,1:num_p)*ones(size(SIG(1:num_p,1:num_p),1),1);
% replace 0s -> machine epsilon
singular_vals(singular_vals<10*eps)=eps;

% Get last columns of V (corresponding to singular vals)

% Case 1: if singular values are eps only:
VT_eps=[];
V_eps=[];
VT_eps=VT(:,find(singular_vals<=10*eps));
V_eps=VT_eps.';

% Case 2: if singular values are > 3 decades below the next highest value
VT_small=[];
V_small=[];
for i=2:length(singular_vals)
    if log10(singular_vals(i))<log10(singular_vals(i-1))-3
        VT_small=VT(:,i:end);
        V_small=VT_small.';
        break
    end
end


D(abs(D)<eps)=nan;

% Sensitivity plot
figure('units','centimeters','position',[0,0,8.78,6.6])
h2=semilogy(t,abs(D),'x--','LineWidth',1.5);
L2=legend(h2,string(p),'Location','northeastoutside');
hold on
L2.AutoUpdate='off';
xlabel('t (s)')
ylabel('|\partial \Delta M / \partial \phi|')
set(gca, 'FontSize', fsz, 'LineWidth', alw)
set(gcf,"Color",'w')
ylim([max(D,[],'all')./10e5 max(D,[],'all').*10])
yticks([0.01 1 100])
xlim([0.9,2.9])

hold off

% Identifiability signature
figure('units','centimeters','position',[0,0,8.78,6.6])
subplot(2,1,1)
pEVs=semilogy(1:1:length(singular_vals),abs(singular_vals),'o','Linewidth',alw);
xticks([1:1:length(singular_vals)])
xtickformat('%i')
ylabel('singular values')
yticks([0.000001 0.01 100])
set(gca, 'FontSize', fsz, 'LineWidth', alw)

V_small(abs(V_small)<0.001)=nan; % remove negligible V values from plot fr clarity

subplot(2,1,2)
if ~isempty(V_small)
    vector_index=length(p)-size(V_small,1);
    for i=1:size(V_small,1)
        vector_index = vector_index +1;
        h(i)=semilogy(abs(V_small(i,:)),'o','Linewidth',alw);
        legend_string(i,:)=['vector ' num2str(vector_index)];
        xlim([1 4])
        hold on
    end
    xticks([1:1:length(p)])
    %xticklabels(string(p))
    xticklabels(["CBF", "ATT", "R1_b", "R1_t"]);
    %set(gca,'TickLabelInterpreter','latex')
    ylabel('V vectors')
    L1=legend(h, legend_string,'AutoUpdate','off');
end
set(gca, 'FontSize', fsz, 'LineWidth', alw)
set(gcf,'Color','w')
hline = refline(0, 0);
hline.Color = 'k';
hline.LineStyle = '--';
L1.AutoUpdate="on";
ylim([0.001 1])

function [f,TA,R1p,R1e,R2p,R2e]=fixparameters(f,TA,R1p,R1e,R2p,R2e)
% Uncomment the fixed parameters:
% R1p  = 1/1.650;
% R1e = 1/1.330;
% f=48*(1/100)*(1/60);
% TA = 1.570;
end