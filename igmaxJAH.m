% JAH 20180426
% Copied from igmax.m, modeling code for the
% DSGE MODEL SMALL OPEN ECONOMY
% CALIBRATED FOR THE KSA ECONOMY
% MODELLING THE ECONOMY WITH DIFFERENT ENERGY SOURCES AND WITH 
% SUBSIDIES FOR DOMESTIC CONSUMPTION AND PUBLIC INVESTMENT IN RENEWABLES
% from the paper "Oil Subsidies and Renewable Energy in Saudi Arabia: A
% General Equilibrium Approach"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Policy Variables
POSUB0 = 0.5;       % posubexc; exogenous value of subsidized oil price
PG0 = 0.6563;       % beta/alpha*posubexg; price of gas
PRNW0 = 1.5625;     %gamma/alpha*posubexg; price of renewable energy
igexg = 0;          %0.00048901051(no subsidies);%0.003327307 (increasing RNW cost);%0.00200347(constant RNW cost);
% above is used to control the % renewable
laborShare = 0.5;   % split of transfers of govt to the 2 social classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
alpha=0.32;         %share of oil in electricity production 
beta=0.42;          %share of gas in electricity production
gamma=1;            %share of renewables in electricity 
a=0.33230;          %share of electricity in energy services
lambda=-0.25786;    %parameter governing the inverse of elasticity of substitution btw electricity and oil in energy services
theta=0.58221;      %labor elasticity
niu=-1.381;         %-1.38;%0.042;%inverse of elasticity of substitution btw capital and energy services in the production function
b=0.00003;          %share of energy services in the production function **************************
d=0.0028;           %weight of household energy services in the utility function
sigmac=-1/3;        %inverse elasticity of substituion between consumption and energy services
sigma=2;            % intertemporal elasticity of substitution
omega=0.0735;       %weight of labor in the utility function
tau=1;              %Frisch elasticity of substitution
phi=0;              %phi=parameter capital adjustment cost
betadisc=0.96;      %betadisc=discount factor
delta=0.1;          %depreciation rate of private capital
nu=0.05;            %depreciation rate of public capital
poaverage=1.6371;   %inconditional average of oil price
A=.0328;            %A=parameter transforming public capital into renewables
rho=0.9;            %rho=persistence parameter in the oil price stochastic process
RBAR=0.04;          %international interest rate
bss=0.16;           %0.117;%trade balance over GDP **************************
psi=0;              %risk premium parameter (0.029 Jesus 0.000742 Uribe)
gexg=0.01398;       %exogenous amount of gas
oexg=0.4819;        %exogenous amount of oil
posubexg=0.5;       %posubexg=exogenous value of subsidized oil price
integcost=0.5;      %degree of the cost of renewable integration into the system
% New Parmaters (JAH)
targRNW = 9.5; % GW target value of renewables for computing govt utility
govUtilWs = [0,0,0]; % need to set these; averaging weights for government utility function

% Initial Values at Time t=0
OE0=0.0305; % volume oil used for electricity production
N0=1; % labor for final goods & services
KP0=5.6246; % capital stock? is this k_t in the paper, i.e. eqn 2
C0=.9363; % consumption of final goods & services
E0=alpha*OE0+beta*gexg+gamma*A*igexg/nu; % electricity produced
G0=gexg; % exogenous amount of gas
RNW0=A*igexg/nu; % renewable energy produced; Rbar_t in the paper
PE0=1/alpha*posubexg; % electricity price
OS0=E0*((1-a)/(alpha*a))^(1/(1-lambda)); % oil used to produce energy services
S0=(a*E0^lambda+(1-a)*OS0^lambda)^(1/lambda); % energy services produced
PS0=PE0/a*(E0^(1-lambda))*S0^(lambda-1); % price of energy services
SF0=KP0*(((RBAR+delta)*b)/(PS0*(1-b)))^(1/(1-niu)); % energy services used to produce final goods & services
W0=theta*N0^(theta-1)*((1-b)*KP0^niu+b*SF0^niu)^((1-theta)/niu); % wages; partial derivative of Y_t wrt nu_t with nu_t=1
R0=RBAR+delta; % related to interest rate?
SH0=S0-SF0; % level of household energy services
RSTAR0=RBAR; % international interest rate
TR0=PRNW0*RNW0+POSUB0*(OS0+OE0)+poaverage*(oexg-OS0-OE0)+PG0*gexg-igexg; % government transfers to households
PO0=poaverage; % inconditional average of oil price (international price)
B0=-bss/RBAR*N0^theta*((1-b)*KP0^niu+b*SF0^niu)^((1-theta)/niu); % level of foreign bonds assetrs
O0=oexg; %exogenous amount of oil
KG0=igexg/nu; % stock of public capital, arising from public investment
Y0=N0^theta*((1-b)*KP0^niu+b*SF0^niu)^((1-theta)/niu); % output of final goods sector

% Create Vectors of Parameters & Initial Values
X0=[OE0,E0,G0,RNW0,PE0,POSUB0,PG0,PRNW0,S0,OS0,PS0,N0,KP0,SF0,W0,R0,SH0,C0,RSTAR0,TR0,PO0,B0,O0,KG0,Y0];
param=[alpha,beta,gamma,a,lambda,theta,niu,b,d,sigmac,sigma,omega,tau,...
    phi,betadisc,delta,nu,poaverage,A,rho,RBAR,bss,psi,gexg,oexg,igexg,posubexg,integcost];

% Minimize to get Steady State of Initial Conditions
 options=optimset('Display','off','MaxFunEvals',100000,'MaxIter',10000);
[XSS,f]=fsolve('ss1',X0,options,param);  % uses Levenberg-Marquadt algo ???

% Extract the Variables
OE = XSS(1);    % amount of oil used for electricity generation
E = XSS(2);     % amount of electricity produced
G =	XSS(3);     % amount of gas used for electricity generation
RNW = XSS(4);   % renewable energy produced
PE = XSS(5);    % electricity price
POSUB = XSS(6); % domesticl oil price
PG = XSS(7);    % price of gas
PRNW = XSS(8);  % price of renewable energy
S = XSS(9);
OS = XSS(10);
PS = XSS(11);   % price of energy services
N = XSS(12);    % wages used to produce final goods
KP = XSS(13);
SF = XSS(14);   % energy services used to produce final goods
W = XSS(15);    % wages used to produce final goods
R = XSS(16);    % return on investing private capital
SH = XSS(17);   % household consumption of energy services
C = XSS(18);    % domestic consumption of final goods
RSTAR = XSS(19);
TR = XSS(20);
PO = XSS(21);   % international oil price
B = XSS(22);
O = XSS(23);
KG = XSS(24);   % stock of private capital
Y = XSS(25);    % output of final goods sector
GDP = Y+PO*(O-OE-OS);
WELFAREINI = 1/(1-sigma)*(C^sigmac+d*SH^sigmac)^((1-sigma)/sigmac);

% Actors' Utilities Calculations
utilElect = PE*E - POSUB*OE - PG*G - PRNW*RNW; % slide 5
% NOPITY NOPE NOPE NOPE elecsUtil = (beta*PO/alpha - PG)*G + (PO/alpha - PRNW)*RNW; % slide 6
utilGoods = Y - W*N - R*KG - PS*SF; % slide 8
utilHouse = ((C^sigmac + d*SH^sigmac)^(1-sigma))/(1-sigma); % slide 9
utilHousL = laborShare*houseUtil; 
utilHousC = (1-laborShare)*houseUtil;
utilGovnt = govUtilWs(1)*WELFAREINI - govUtilWs(2)*abs(RNW-targRNW) - govUtilWs(3)*(poaverage/PO);
