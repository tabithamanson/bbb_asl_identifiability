% ALL UNITS ARE IN SECONDS NOW

M0=5400*0.85;%4590*0.85; %average from 23 control participants DPRCX
%M0=100.*0.85; % to get everything as a percentage of M0a
vp=0.05; %https://www.sciencedirect.com/topics/nursing-and-health-professions/brain-blood-volume
R1p=1/1.650; %1.86 for zhao et al 2007 calculated value for HCT 0.36 Y=0.77 (capillaries)
R1e=1/1.330;
f=48*(1/100)*(1/60); %ml / s / g
F=f;
TA=1.57;%1.570; %s
kw=140*(1/60); % /s
Texch=1/kw; %s
dA=TA;%1.57
d=dA+Texch;
R2p=1/0.110;
R2e=1/0.077;
lambda=0.9;

% for biexponential model:
T2p=0.110;%106 for zhao et al calculated value
T2e=0.077;

%Sp=0.45*0.02*(1/100)*(M0/0.85); %TI=2100
%Se=0.55*0.02*(1/100)*(M0/0.85); %TI=2100
%  Se=5.5; %TI=2100
%  Sp=2.1; %TI=2100
% Se=3.5; %TI=3100
% Sp=0.11; %TI=3100
%C=0;

%Se=9; %2100
%Sp=3.7;

 Se=5.9; %3100
 Sp=0.2;

Scsf=(Se+Sp)/100;
T2csf=1.7;
Acsf=(Se+Sp)/1000;
kcsf=1/60;
%vcsf=0.1;


