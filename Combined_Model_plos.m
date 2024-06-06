clear
clc
tic
%% %%%%%%%%%%%% (c)2021, Ahmed Aldohbeyb & Kvein Lear (10/7/2021) 
%%%%%%%%%%% 
% A single-compartment HH-type model that combining the effect of cooperative
% Na+ channels and dynamical reversal potential, which we called the 
% Combined Model. The model is based on Pospischil et al (2008) paper, which 
% described different type of neurons and their VGICs. The dynamical reversal potential was 
% adopted from Cressman et al (2009I) full model and cooperative Na+ channels 
% were used as described in Huang et al (2012) model. 
% Note that some of the parameters was based on Cressman et al(2009) published 
% model in ModelDB where several corrections was made in order to accurately 
% reproduce their results. 
%% cell type as described in Pospischil et al (2008): 
% 1) regular-spiking (RS) excitatory (pyramidal) neuron & inhibitory neuron
% 2) fast-spiking (FS) neuron,
% Note: changing some the parameters such as g, VT, or tau can reproduce the response of different cells that fall
% under the same categorize. for example, depending on the choise of Na, Kd, Ks, and leak channels
% parameters, you can get the response from a RS excitatory or RS inhibitory cell. 
% So, this code can reproduce some of the figures in Pospischil et al
% (2008) paper by choosing the cell type and setting dyn_con = 0 & KJ or p = 0.
%% neuron types 
% 1: RS (exc (PC)) based on somatosensory cortex in vitro fig.2a in Pospischil et al (2008),
% 2: RS(inh), based on somatosensory cortex in vitro fig.2b in Pospischil et al (2008),
% 3: FS (inh), based on ferret visual cortex in vitro Fig.3 in Pospischil et al (2008), 
% 4: FS (inh), based on somatosensory cortex in vitro Fig.4 in Pospischil et al (2008), 
% model parameter for each cell type. Each coulum represent a neuron type and
% the rows represent the model parameters as following: 1)Na+ conductance,
% 2)g_Kd conductance, 3) g_Ks conductance, 4) g_leak conductance 5) E_leak, 
% 6) VT 7)time const. for g_Ks, 8) model dimension in 
% um (L = d). Units: conductance (mS/cm2), potential (mV), time (ms), and
% length and depth in um.
% dyn_cond = 1 : dynamical reversal potential & 0 : steady-state value (fixed) 
tic
dyn_cond = 1; 
type = 3; 
para = [56 10 50 58;6 2.1 10 3.9;7.5E-2 9.8E-2 0 7.87E-2;...
 2.05E-2 1.33E-2 1.5E-1 3.8E-2;-70.3 -65 -70 -70.4; ...
 -56.2 -67.9 -63 -57.9;608 934 4000 502;61.4 61.8 67 56.9]; 
KJ = 0;
p = 0.0;                                                                   % fraction of coop. chs' (%)
KJ = 0;                                                                    % coupling strength (mV);
%% Parameters
dt = 1E-3;                                                                 % time step (ms)
T_max = 2.0E3;                                                             % max time (ms) Note: increase the number for multiple spike trains
t = 0:dt:T_max;                                                            % time vector (ms) 
N = length(t);                                                             % number of steps 
A = (para(8,type)*1e-4)^2;                                                 % Area para(10,type) is in um, so x1E-4 to have  A in cm^2
C = A*1;                                                                   % Membrane capacitance (1 uF/cm^2)
gNa = para(1,type);                                                        % Na+ max conductance (mS/cm^2)
gKd = para(2,type);                                                        % K+ (delayed-rectifier) max conductance (mS/cm^2) 
gKs = para(3,type);                                                        % K+ (Slow) max conductance (mS/cm^2)
gL = para(4,type);                                                         % Leak max conductance (mS/cm^2) 
EL = para(5,type);                                                         % Leak reversal potential (mV) 
VT = para(6,type);                                                         % Variable to adjusts spike threshold
tau_nsMax = para(7,type);                                                  % maximum time const. for slow K+ ch's
deg = 34.1;                                                                % Tempreture °C
V_thermal = 1000*(deg + 273)*(1.38E-23/1.6E-19);                           % Thermal voltage 
B = 7;                                                                     % Ratio of intracellular to extracellular volume of the cell
Ps = 1.125 ;                                                               % Pump strength (mM/ms)
G_glia = 200/3;                                                            % Strength of glial uptake (mM/ms)
epsilonK =4.0/3.0;                                                         % Diffusion constant (1/ms)
kbath = 4.0;                                                               % Steady state extracellular potassium concentration (mM)
gamma = 0.044495;                                                          % factor to convert current density to d[conc]/dt (mM.cm2/?coul)
tau =1000;                                                                 % to convert mM/s to mM/ms
Na_i = zeros(1,N);Na_o=zeros(1,N);                                         % preallocate vectors for ion concentration 
K_i=zeros(1,N);K_o= zeros(1,N);
Na_i(1) = 18.0;                                                            % Initial intracellular Na+ conc at the channels (mM)
% Na_p = Na_i;                                                             % Intial intracellular Na+ conc at the pump (mM)
Na_o(1) = 145;                                                             % Initial extracellular Na+ conc (mM)
K_i(1) = 140 ;                                                             % Initial intracellular K+ conc (mM)
K_o(1) = 3.7;                                                              % Initial extracellular K+ conc (mM)
kna = 87.5;                                                                % Na+ dissociation constant (mM)
kca = 1.38;                                                                % Ca2+ dissociation constant(mM)
ksat = 0.1;                                                                % the saturation factor
GM = 0.35;                                                                 % voltage dependence factor
q10 = 3^((deg - 37)/10); 
% D_Na = 3e-6;
% delta_X = 1E-4;
% epsilonNa = 2*D_Na./(delta_X)^2;
ENa = zeros(1,N);                                                          % preallocate vectors for reversal potential 
EK = zeros(1,N);
ENa(1) = V_thermal*log(Na_o(1)./Na_i(1));                                  % Initial Na+ reversal potential (mV)
EK(1) = V_thermal*log(K_o(1)./K_i(1));                                     % Initial K+ reversal potential Na+ (mV) 
%% Intial Conditions
V = zeros(1,N);                                                            % preallocate vectors for Vm and gating variables
Vcoop=V;m=V;
mc=V;h=V;
n=V;ns=V; 
V(1) = EL;                                                                 % initial membrane potential
m(1) = alphaM(V(1),VT)/(alphaM(V(1),VT) + betaM(V(1),VT));                 % independent Na+ channels activation 
mc(1) = m(1);                                                              % cooperative Na+ channels activation 
h(1) = alphaH(V(1),VT)/(alphaH(V(1),VT) + betaH(V(1),VT));                 % Na+ channels inactivation 
n(1) = alphaN(V(1),VT)/(alphaN(V(1),VT) + betaN(V(1),VT));                 % K+ (delayed-rectifier) channels activation 
ns(1) = Ns_inf(V(1));                                                      % K+ (slow) channels activation 
INai = zeros(1,N);INac = INai;INa = INai;                                  % preallocate vectors for ch's current 
IKd = zeros(1,N);IKs = IKd;IK = IKd;
ILeak = zeros(1,N);Ipump=ILeak;Iglia=ILeak;IKdiff =ILeak;
%% current pulse & synaptic parameters
I = 0;                                                                     % Vector for current pulse
I_app = 2.5;
PL_T = round(1000/dt);                                                     % length of the current pulse 
T_bw_pulses = round(5000/dt);                                              % Time between current pulses 
I_on = round(500/dt);                                                     % on time for 1st current pulse
I_off = I_on + PL_T;                                                       % off time for 1st current pulse
I_on2 = I_off + T_bw_pulses;                                               % on time for 2nd current pulse
I_off2 = I_on2 +PL_T; % off time for 2nd current pulse
I_on3 = I_off2 + T_bw_pulses; % on time for 3rd current pulse
I_off3 = I_on3 +PL_T; % off time for 3rd current pulse
I_on4 = I_off3 + T_bw_pulses; % on time for 4th current pulse
I_off4 = I_on4 + PL_T; % off time for 4th current pulse
I_on5 = I_off4 + T_bw_pulses; % on time for 4th current pulse
I_off5 = I_on5 + PL_T; % off time for 4th current pulse
%% the main loop
for i = 1:N-1
%% Current through each voltage-gated ion channels & Membrane potential 
INai(i) = (1-p)*gNa*m(i)*h(i)*(ENa(i) - V(i));                             % Independent Na+ channels current
INac(i) = p*gNa*mc(i)*h(i)*(ENa(i) - V(i));                                % Cooperative Na+ channels current
INa(i) = INai(i) +INac(i);                                                 % total Na+ current
IKd(i) = gKd*n(i)^4*(EK(i) - V(i));                                        % delayed-rectifier K+ current
IKs(i) = gKs*ns(i)*(EK(i) - V(i));                                         % slow non-inactivating K+ current 
IK(i) = IKd(i) +IKs(i);                                                    % total K+ current 
ILeak(i) = gL*(EL-V(i));                                                   % leak channel current
% membrane potential (4RK method)
Vrk(1) = (dt/C)*A*(INa(i) + IK(i) + ILeak(i) + I); 
kV = V(i) +Vrk(1)/2;
Vrk(2) = 0.5*(dt/C)*A*((1-p)*gNa*m(i)*h(i)*(ENa(i) - kV) + p*gNa*mc(i)*h(i)*(ENa(i) -kV)...
 + gKd*n(i)^4*(EK(i) - kV)+ gKs*ns(i)*(EK(i) - kV) + gL*(EL-kV) + I); 
kV = V(i) +Vrk(2)/2;
Vrk(3) = 0.5*(dt/C)*A*((1-p)*gNa*m(i)*h(i)*(ENa(i) - kV) + p*gNa*mc(i)*h(i)*(ENa(i) -kV)...
 + gKd*n(i)^4*(EK(i) - kV)+ gKs*ns(i)*(EK(i) - kV) + gL*(EL-kV) + I); 
kV = V(i) +Vrk(3);
Vrk(4) = (dt/C)*A*((1-p)*gNa*m(i)*h(i)*(ENa(i) - kV) + p*gNa*mc(i)*h(i)*(ENa(i) - kV)...
 + gKd*n(i)^4*(EK(i) - kV)+ gKs*ns(i)*(EK(i) - kV) + gL*(EL-kV) + I); 
V(i+1) = V(i) + (Vrk(1) + 2*Vrk(2) + 2*Vrk(3) + Vrk(4))/6;
%% Conc. & Reversal potential (Dynamical (dyn_con = 1) or Fixed (dyn_con = 0))
if(dyn_cond ==0)                                                           % if fixed reversal potential use conc. SS value
ENa(i+1) = V_thermal*log(Na_o(1)./Na_i(1));
EK(i+1) = V_thermal*log(K_o(1)./K_i(1));
else                                                                       % if dynamical reversal potential 
ENa(i+1) = V_thermal*log(Na_o(i)./Na_i(i));                                % Na+ reversal potential (mV)
EK(i+1) = V_thermal*log(K_o(i)/K_i(i));                                    % K+ reversal potential (mV)
Ipump(i) = (Ps/(1.0+exp((25.0-Na_i(i))/3.0)))*(1/(1+exp(5.5-K_o(i))));     % Na+/K+ pump (mM/s)
Iglia(i) = G_glia/(1.0+exp((18.0-K_o(i))/2.5));                            % glial capacity to remove excess K+ from the extracellular space (mM/s)
IKdiff(i) = epsilonK*(K_o(i)-kbath);                                       % diffusion of potassium away from the local extracellular micro-environment (mM/s)
% INadiff(i) = epsilonNa*(Na_i(i)-Na_p(i));                                % diffusion of Na from the channels to the Na/K pump
Na_o(i+1) = 145 - B*(Na_i(i)-Na_i(1));                                     % Extracellular Na+ concentration (mM)
K_i(i+1) = 140 +(Na_i(1) - Na_i(i));                                       % Intracellular K+ concentration (mM)
% 4RK for [Na+]i & [K+]o
% Intial values 
INa0 = INa(i);
Iglia0 = Iglia(i);
IKdiff0 = IKdiff(i);
% INadiff0 = INadiff(i);
IK0 = IK(i);
INaK = Ipump(i);
dt0=dt;
for j = 1:4                                                                % Calculate the 4 RK coefficients for each ion species
k_K(j)= dt0*(-gamma*B*IK0 - 2*B*INaK -Iglia0-IKdiff0)/tau; 
k_Na(j) = dt0*(gamma*INa0-3*INaK)/tau;
% k_Nap(j) = dt0*(-3*INaK + INadiff0)/tau;
if(j==3)                                                                   % for the fourth coefficients, the values are multiplied by 1
Na_i0 = Na_i(i) + k_Na(j);
% Na_p0 = Na_p(i) + k_Nap(j);
K_o0 = K_o(i) + k_K(j);
dt0 = dt; 
else 
Na_i0 = Na_i(i) + k_Na(j)/2;
% Na_p0 = Na_p(i) + k_Nap(j)/2;
K_o0 = K_o(i) + k_K(j)/2;
dt0 = dt/2;
end
ENa0 = V_thermal*log(Na_o(i)./Na_i0); 
INa0 = (1-p)*gNa*m(i)*h(i)*(ENa0 - V(i)) + p*gNa*mc(i)*h(i)*(ENa0 - V(i)); 
EK0 = V_thermal*log(K_o0./K_i(i)); 
Iglia0 = G_glia/(1.0+exp((18.0-K_o0)/2.5)); 
IKdiff0 = epsilonK*(K_o0-kbath); 
% INadiff(i) = epsilonNa*(Na_i0-Na_p0);
IK0 = gKd*n(i)^4*(EK0 - V(i)) + gKs*ns(i)*(EK0 - V(i));
INaK = (Ps/(1.0+exp((25.0-Na_i0)/3.0)))*(1/(1+exp(5.5-K_o0))); 
end
K_o(i+1)=K_o(i) + (k_K(1) + 2*k_K(2) +2*k_K(3) +k_K(4))/6;                 % Extracellular K+ concentration (mM)
Na_i(i+1)=Na_i(i) + (k_Na(1) + 2*k_Na(2) +2*k_Na(3) +k_Na(4))/6;           % Intracellular Na+ concentration (mM)
% Na_p(i+1)=Na_p(i) + (k_Nap(1) + 2*k_Nap(2) +2*k_Nap(3) +k_Nap(4))/6;     % Na+ concentration at the pump(mM)
end
%% gating variables (4RK method for differential eq's)
Vcoop(i) = V(i) + KJ*mc(i)*h(i);                                           % V for Na+ coopetive channels
%% Na+ gating varibles 
mc_inf = alphaM(Vcoop(i),VT)/(betaM(Vcoop(i),VT) + alphaM(Vcoop(i),VT)); 
tau_mc = 1/(alphaM(Vcoop(i),VT)+betaM(Vcoop(i),VT));
mc(i+1) = mc_inf^3 + (mc(i)-mc_inf^3)*exp(-dt/tau_mc);
m_inf = alphaM(V(i),VT)/(betaM(V(i),VT) + alphaM(V(i),VT));
tau_m = 1/(alphaM(V(i),VT)+betaM(V(i),VT));
m(i+1) = m_inf^3 + (m(i)-m_inf^3)*exp(-dt/tau_m);
fh = @(h0,V,dt) dt*(0.128*exp(-(V-VT-17)/18)*(1-h0) - h0*4./(exp(-(V-VT-40)/5)+1));
h(i+1) = RK4_function(fh,dt,h(i),V(i));
%% K+ gating varibles
fn = @(n0,V,dt) dt*(-0.032*(V-VT-15)/(exp(-(V-VT-15)/5) - 1)*(1-n0) - n0*0.5*exp(-(V-VT-10)/40));
n(i+1) = RK4_function(fn,dt,n(i),V(i));
fns = @(ns0,V,dt) dt*((1/(1+exp(-(V+35)/10))) - ns0)/(tau_nsMax/(3.3*(exp((V+35)/20)) + exp(-(V+35)/20)));
ns(i+1) = RK4_function(fns,dt,ns(i),V(i));
%% current pulses
if (i>I_on)&&(i<I_off)                                                    % 1st pulse
 I = I_app; 
 elseif(i>I_on2)&&(i<I_off2) % 2nd pulse
 I =1.5*I_app;
 elseif(i>I_on3)&&(i<I_off3) % 3rd pulse
 I = 2*I_app;
 elseif(i>=I_on4)&&(i<I_off4) % 4th pulse
 I = 2.5*I_app; 
 elseif(i>=I_on5)&&(i<I_off5) % 4th pulse
 I = 3*I_app;
 else 
 I = 0;
 end
%% check if V(i) becomes a complex number & break if it is a complex number or V> or < +- 300 mV
 tf = isreal(V(i)); 
if(tf == 0) || (abs(V(i))>300) 
mass = [' Error occured. V = ',num2str(round(V(i))),' mV at t = ',num2str(t(i)),' ms'];
disp(mass)
 break
end
end
figure(1)
plot(t./1000, V, 'LineWidth',1.5)
ylabel('V$_m$ (mV)','FontSize',16,'Interpreter','latex')
xlabel('Time (sec)','FontSize',16,'Interpreter','latex')
title('Membrane Potnetial','FontSize',16,'Interpreter','latex')
figure(2)
plot(t./1000, ENa, 'LineWidth',1.5)
ylabel('E$_{Na}^+$ (mV)','FontSize',16,'Interpreter','latex')
xlabel('Time (sec)','FontSize',16,'Interpreter','latex')
title('Reversal Potential (Na$^+$)','FontSize',16,'Interpreter','latex')
figure(3)
plot(t./1000, EK, 'LineWidth',1.5)
ylabel('E$_{K}^+$ (mV)','FontSize',16,'Interpreter','latex')
xlabel('Time (sec)','FontSize',16,'Interpreter','latex')
title('Reversal Potential (K$^+$)','FontSize',16,'Interpreter','latex')
time = toc/60; % Simulation time in minutes 
%% rate functions
% alpha
function aM = alphaM(V,VT)
aM = -0.32*(V-VT-13)/(exp(-(V-VT-13)/4) - 1);
end
function aH = alphaH(V,VT)
aH = 0.128*exp(-(V-VT-17)/18);
end
function aN = alphaN(V,VT)
aN = -0.032*(V-VT-15)/(exp(-(V-VT-15)/5) - 1);
end
% beta
function bM = betaM(V,VT)
bM = 0.28*(V-VT-40)/(exp((V-VT-40)/5) - 1);
end 
function bH = betaH(V,VT)
bH = 4./(exp(-(V-VT-40)/5)+1);
end
function bN = betaN(V,VT)
bN = 0.5*exp(-(V-VT-10)/40);
end
function Ns_i = Ns_inf(V)
Ns_i= 1/(1+exp(-(V+35)/10));
end
function tns_i = tau_ns(V,tau_nsMax)
tns_i= tau_nsMax/(3.3*(exp((V+35)/20)) + exp(-(V+35)/20));
end
%% 4th RK for channels variables
function RK4g= RK4_function(f,dt,a0,V)
dt0 =dt;
m1 = a0;
for j = 1:3
 k(j) = f(m1,V,dt0);
 m1 = a0 + k(j)/2;
 dt0 = dt/2;
end
m1 = a0 + k(j);
 dt0 = dt;
 k(4) = f(m1,V,dt0);
RK4g = a0 + (k(1) + 2*k(2) +2*k(3) +k(4))/6;
end
