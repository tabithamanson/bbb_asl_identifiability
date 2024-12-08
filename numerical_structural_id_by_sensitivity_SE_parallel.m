% Structural identifiability analysis using the sensitivity matrix. Parallel 2CXM, single-echo scan.
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
% clear workspace before running 

clear workspace
sympref('TypesetOutput',true);
syms y(t) x(t)
syms M0 f F TA R1p kw tau R1e vp t R2p R2e k positive real rational
syms s positive

alw = 2;    % AxesLineWidth
fsz = 12;      % Fontsize

%make a struct of free parameters
p=[f TA R1p R1e kw];% M0 tau];
m_max=4; % What is this?
% 5 for parallel_2CXM, 
% 7 for series_2CXM or 
% 9 for buxton heaviside, 4 for c=100

% initial condition
 x0=[0; 0];

vars=[p];

num_p=length(p);
num_x=size(x0,1);

% scale factor big "C"
C=[1 1];

y(t,p)=C*x;

%SE t_array
t_array=[500:400:2900]./1000; %s

num_t=length(t_array);

tau=0.400;
run('literature_vals.m')
c=100;
heaviside_approx=1;
syms y_subs(t)
y_subs(t)=subs(y);   
p_subs=subs(p);

% create matrix of central differences about current literature values:

%left parameter values = 0.99*param vals
%right param vals = 1.01*param vals
change(1)=+Inf; % first change in S entries (from "no value" to something)
change_max(1)=+Inf;
m=0;
p_del=+Inf;
while (sum(p_del<10*eps) < 1 && m<m_max)
m=m+1;
del_L=-(1e-2)*1*10^(-(m-1));
del_R=-del_L;
    for k=1:(num_p) %num parameters
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
[f,TA,R1p,R1e,R2p,R2e,kw]=fixparameters(f,TA,R1p,R1e,R2p,R2e,kw);

[Mt_L,Mp_L,Me_L] = my_solve_parallel_2CXM_numerical(c,TA,tau,R1p,R1e,f,kw,M0,t_array,heaviside_approx);


     t=t_array;

    M(:,1,k)=Mt_L;

    run('literature_vals.m')
    end
    
        for k=1:(num_p) %num parameters
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
     [f,TA,R1p,R1e,R2p,R2e,kw]=fixparameters(f,TA,R1p,R1e,R2p,R2e,kw);
     [Mt_R,Mp_R,Me_R] = my_solve_parallel_2CXM_numerical(c,TA,tau,R1p,R1e,f,kw,M0,t_array,heaviside_approx);
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
   
   % remove small values of v

   % for j=1:size(V_small,1)
   %     [V_sorted, I_v]=sort(abs(V_small(j,:)),'descend');
   %     for i=2:size(V_small,2)         
   %          if log10(V_sorted(i))<log10(V_sorted(i-1))-3
   %              V_small(j,I_v(i:end))=0;
   %              break
   %          end
   %     end
   % end

    V_small(abs(V_small)<0.001)=0;

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

%get last columns of V (corresponding to singular vals) 

% case1: if singular values are < eps:
VT_eps=[];
V_eps=[];
VT_eps=VT(:,find(singular_vals<=10*eps));
V_eps=VT_eps.';

% %Case 2: if singular values are > 3 decades below the next highest value
VT_small=[];
V_small=[];
    for i=2:length(singular_vals)
        if log10(singular_vals(i))<log10(singular_vals(i-1))-3
            VT_small=VT(:,i:end);
            V_small=VT_small.';
            break
        end
    end
V_small(abs(V_small)<0.001)=0;

%plot identifiability signature (case 2)
% % figure(2)
% % subplot(2,1,1)
% % pEVs=plot(1:1:length(singular_vals),log10(singular_vals),'o','Linewidth',alw);
% % title(["Identifiability Signature - all v"])
% % xticks([1:1:length(singular_vals)])
% % xtickformat('%i')
% % ylabel('log(eigenvalues)')
% % % hline = refline(0, 0);
% % % hline.Color = 'k';
% % set(gca, 'FontSize', fsz, 'LineWidth', alw)
% % set(gcf,'Color','w')
% % hold on
% % pEPS=plot(1:1:length(singular_vals),ones(1,length(singular_vals)).*log10(eps),'--');
% % legend([pEPS],'log(machine-zero)');
% % 
% % subplot(2,1,2)
% % if ~isempty(V_eps)
% % vector_index=length(p)-size(V_eps,1);
% % for i=1:size(V_eps,1)
% % vector_index = vector_index +1;
% % h(i)=plot(V_eps(i,:),'o','Linewidth',alw);
% % legend_string(i,:)=['vector ' num2str(vector_index)];
% % hold on
% % %legend(h, legend_string)%,'AutoUpdate','off')
% % end
% % end
% % xticks([1:1:length(p)])
% % xticklabels(string(p))
% % ylabel('V vectors')
% % set(gca, 'FontSize', fsz, 'LineWidth', alw)
% % set(gcf,"Color",'w')
% % hline = refline(0, 0);
% % hline.Color = 'k';
% % hline.LineStyle = '--';
% % 


%hold off
D(abs(D)<eps)=eps;

figure('units','centimeters','position',[0,0,8.78,6.6])
h2=semilogy(t,abs(D),'x--','LineWidth',1.5);
L2=legend(h2,string(p),'Location','northeastoutside');
hold on
L2.AutoUpdate='off';
%xline(gca,0.5:0.400:2.900,'k--')
%xlim([0.4 3.2])
% pEPS=plot(t,ones(1,length(t)).*(eps),'--');
% legend_string2=[string(p),"machine-zero"];
% 
% legend([h2;pEPS],legend_string2)
xlabel('t (s)')
ylabel('|\partial \Delta M / \partial \phi|')
set(gca, 'FontSize', fsz, 'LineWidth', alw)
set(gcf,"Color",'w')
ylim([1*10^(-6) max(D,[],'all').*10])
xlim([0.9,2.9])
yticks([0.0001 0.01 1 100])
hold off


figure('units','centimeters','position',[0,0,8.78,6.6])
%set(gcf, 'Position',  [100, 100, 500, 400])
subplot(2,1,1)
pEVs=semilogy(1:1:length(singular_vals),abs(singular_vals),'o','Linewidth',alw);
%title(["Identifiability Signature - small V removed"])
xticks([1:1:length(singular_vals)])
xtickformat('%i')
ylabel('singular values')
yticks([0.000001 0.01 100])

set(gca, 'FontSize', fsz, 'LineWidth', alw)
% set(gcf,'Color','w')
% hline = refline(0, 0);
% hline.Color = 'k';

% hold on
% pEPS=plot(1:1:length(singular_vals),ones(1,length(singular_vals)).*log10(eps),'--');
% legend([pEPS],'log(machine-zero)')
%legend
V_small(abs(V_small)<0.001)=nan;
subplot(2,1,2)
if ~isempty(V_small)
vector_index=length(p)-size(V_small,1);
for i=1:size(V_small,1)
vector_index = vector_index +1;
h(i)=semilogy(abs(V_small(i,:)),'o','Linewidth',alw);
legend_string(i,:)=['vector ' num2str(vector_index)];
hold on
end
xticks([1:1:length(p)])
xlim([1 5])
%xticklabels(string(p))
xticklabels(["CBF", "ATT", "R1_b", "R1_t", "k_w"]);
%set(gca,'TickLabelInterpreter','latex')
ylabel('V vectors')
L1=legend(h, legend_string,'AutoUpdate','off');
end
set(gca, 'FontSize', fsz, 'LineWidth', alw)
set(gcf,'Color','w')
xticklabels(["CBF", "ATT", "R1_b", "R1_t", "k_w"]);
hline = refline(0, 0);
hline.Color = 'k';
hline.LineStyle = '--';
L1.AutoUpdate="on";
ylim([0.001 1]);



%For 1 : m columns
%calculate dY(..m)/dTheta(...p+n)

%%
function [f,TA,R1p,R1e,R2p,R2e,kw]=fixparameters(f,TA,R1p,R1e,R2p,R2e,kw)
%[M0 tau F TA R1p R1e R2p R2e] = [1 1000 48*(1/100)*(1/60)*(1/1000) 1570) 1/1650 1/1330 1/110 1/77]
% R1p  = 1/1.650;
% R1e = 1/1.330;
% f=48*(1/100)*(1/60)*(1/1000);
% R2p=1/0.110;
% R2e=1/0.077;
% TA = 1.570;
% kw=140/60;
end