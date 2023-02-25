clear all
close all
clc
% parpool(16);
%% Folder Initialization

disp('Adding Initialization Folder Files, Please Connect External HDD ....')
home_fold = pwd;
if ismac
    slash = '/';
else
    slash = '\';
end

disp('Adding Machine Metamodel ....')
addpath([home_fold slash 'Machine Metamodel 2']);
disp('Done')

%% Structure Initializtion
% Simulation paramters ----------------------------------------------------
D.N.Nppsw = 140;                     % number of points per switching cycle

% capacitor parameters ----------------------------------------------------
% Maximum DC bus
D.C.Vcdsmax = 1500;

% Conductivity
D.C.Gvbase = 1200;

D.C.aga = 1.1948;
D.C.agb = 0.3619;
D.C.agc = 532.8210;
D.C.agd = 3.9168;

D.C.bga = 0.04067;
D.C.bgb = 1.8335e5;
D.C.bgc = 0.7285;
D.C.bgd = 1.0168e-4;

% Maximum temperature for 100,000 Hours
D.C.Tmax = 70;
D.C.Tamb = 25;
D.C.DTempmax = D.C.Tmax - D.C.Tamb ;

% ESR
D.C.ESR.beta = 0.002177744479519;
D.C.ESR.alpha1 = 2.565678718044898e-09;
D.C.ESR.alpha2 = 6.521586111532122e-08;
D.C.ESR.gamma1 = 9.695833826880358;
D.C.ESR.gamma2 = 0.255739857707302;
D.C.ESR.vb1 = 5.842874497659989e+02;
D.C.ESR.vb2 = 7.634787419819759e+02;

% Mass
D.C.M.beta = 0.017950077893489;
D.C.M.alpha = 27.551645401982910;
D.C.M.gamma = 0.642829822786177;
D.C.M.vb = 1.785275430393037;

D.C.vrip   = 0.05; % how many times of Dc voltage should be ripple voltage max - allow 5 % ripple

% Heat Sink ---------------------------------------------------------------
% fit equation a/(x)^n1 + b/(x)^n2
% D.HS.a = 0.004867; D.HS.n1 = 1.679;
% D.HS.b = 0.06314; D.HS.n2 = 0.3263;
D.HS.a = 6.33e-3; D.HS.n1 = 1.63;
D.HS.b = 56.28e-3; D.HS.n2 = 0.36;
% fan power
D.HS.f.a = 0.5645; D.HS.f.n1 = 1.209;
D.HS.f.b = 0.2545; D.HS.f.n2 = -0.1992;

D.HS.Thsmax = 100; D.HS.Tjmax = 100;
D.HS.Tamb = 25; D.HS.Rmin = 0.1; D.HS.Rlim = 3; D.HS.Mmax = 1;
% mfit = D.HS.a./(Rthsa).^D.HS.n1 + D.HS.b./(Rthsa).^D.HS.n2;

% % Capacitor Life Time Pararameters
% % Ambient temperature
% D.C.tamb = 20; % ambient temp in degree C
% D.C.dt = 5 ;  % change in temperature
% D.C.L1 = 12000; % life time of capacitor for 60 degree C
% D.C.tmax = 85; % temp in degree C - rateed temperature

% UI core parameters ------------------------------------------------------
% Ferrite core dc inductor sizing parameters
D.L.cM = 5.9851;
D.L.bM = [ 0 100.05 100.05 3.4677e6 8.2537e6 7.3079e7 1.0430e8 ];
D.L.nM = [ 0.24700 0.24673 -1.3215 -1.2423 2.4809 -2.0633 1.5530 ];
D.L.cP = 1.1021e-5;
D.L.bP = [ 0 1.1658e3 1.1658e3 1.1659e3 5.1412e4 4.4344e5 1.2330e6 ];
D.L.nP = [ 0.54482 0.25254 0.17114 0.24906 -0.15241 0.52755 -0.59614 ];

% Switch Parameters -------------------------------------------------------
% Eonn Eoffn Ern are normalized values
% data where Eon and Eoff is not available but Etotal is available, Eon and
% Eoff are assumed equal

% Voltage   Current   Volume        Eoffn        Eonn         Ern
% rating    rating      
sw_data = [ ...
  250        400     1.0644e-04   3.7133e-07   3.7133e-07   7.5e-7    ;
  600        600     5.6916e-05   1.5690e-07   3.7720e-08  1.2781e-07 ;
 1200        600     1.0044e-04   1.8530e-07   1.8845e-07  1.2847e-07 ];

% Conduction Loss Data
% Vsw        Rsw         Vfd         Rfd
sw_data_1 = [ ...
0.7705     0.001174    0.6911     0.002598 ;
0.6000     0.001528    0.5600     0.001700 ;
0.8000     0.002700    0.7300     0.003200 ];


% Thermal Data
% Rthswj-c    Rthdj-c    Rthc-h 
sw_data_2 = [ ...
   0.09         0.20       0.1    ;
   0.065        0.11       0.1    ; % Rthc-h was not available, assumed a typical number
   0.066        0.12       0.02   ];

% Mass Data
% Mass_Package  Gate_Driver
% devide by half because of 2 devices in one package
sw_data_3 = [ ...
   400/2e3        50/4e3 ;
   300/2e3        50/4e3 ;
   400/2e3        50/4e3 ];


D.sw_data = [sw_data sw_data_1 sw_data_2 sw_data_3];

% New Conduction loss data
% V = beta + alpha*I^gamma
% [beta alpha gamma]
% [Switch Diode]
D.swdata2 = [0.7216 0.0022 0.9098 0.4224 0.0493 0.5335; ...
    0.2547 0.0412 0.5306 0.4766 0.0288 0.5411; ...
    0.7363 0.0089 0.8128 0.4857 0.0467 0.5730];

% Leg Current Limits/ Current Density Limits ------------------------------
D.irip  = 1.07; % 7 percent of the fundamental component

D.Jpeak = 10e6; % peak current density of inductor
D.Jpeak_m = 10e6;% peak current density of machine

D.fswmax = 50e3; % removed pk freq limit to 40 kHz
D.fundlim = 200; % 200 Hz fundamental comp



% DC Bus Parameters -------------------------------------------------------
D.Pdc = 1e6; % 200 kW power supply
D.Vdc  = 5e3 ; % fixed
D.vbusrip = 0.001;

% Load Metamodel Paramters ------------------------------------------------
load('SMPMG_fr_asr.mat');
D.SMPMG = SMPMG;
D.tipspmax = 200;  % maximum tip speed
D.SMPMG.l_fac = 1.1; D.SMPMG.d_fac = 1.1; 
D.SMPMG.Ns1max = 1e50; % maximum number of fundamental comp of turns in machine

% Mass Factors
% % for machine
% D.SMPMG.m_fac = 1.5; % EM /Real Mass
% % for submodule 
% D.MfacSM = 1.15;
% % for inductor
% D.MfacI = 1.05;

% for machine
D.SMPMG.m_fac = 1; % EM /Real Mass
% for submodule 
D.MfacSM = 1;
% for inductor
D.MfacI = 1;


% Put a maximum limit of losses -------------------------------------------
D.Pfc = 1.20; % maximum 20% losses

% Maximum allowed mass
D.Mmax = 1e4;

% Modulation Type
% 1 means upper and lower legs use 180 degree shifted triangle wave
% otherwise it means same triangle is uses for the whole mmc simulation
% put it as a gene, why not?
D.ModType = 0;


%% genetic algorithm parameters -------------------------------------------
ngen=3000;
npop=3500;
D.UnitMetric = [1 1];
GAP=gapdefault(numel(D.UnitMetric),0,npop,ngen);
% GAP.dv_act=0;
GAP.ev_pp = 1;
GAP.ev_npg = 32;
GAP.rp_lvl=1;
GAP.rp_gbr = 10;
%%
disp(['Pdc is ' num2str(D.Pdc/1e3) 'kW'])
disp(['Vdc is ' num2str(D.Vdc/1e3) 'kV'])
%% gene description--------------------------------------------------------
%
%         min   max chrm  chrm      par
%         val   val type   id         #     description
GAP.gd= [  
%----------------------- Machine Parameters ------------------------------%

          1e2   8e3   3    1;  ...  %  1 Generator Speed in RPM
           1    20    3    1;  ...  %  2 Torque Density of Machine
           1   1e5    3    1;  ...  %  3 Estimated MMC Loss
          0.5   1     3    1;  ...  %  4 Generator Ns1a
          
%----------------------- MMC Parameters ----------------------------------%

         1e-4  5e-3   3    1;  ...  %  5 Capacitance of each cap in Submodule
           4     20   3    1;  ...  %  6 Number of Submodules in each leg
          50    1e3   3    1;  ...  %  7 Number of cycles per fundamental cycle, increased the minimum required cycles to make te algo accurate
         1e-9  1e-2   3    1;  ...  %  8 Self Inductance
         1e-4  5e-3   3    1;  ...  %  9 Bus Capacitance
         
%----------------------- EE Core Parameters ------------------------------%
         
        5e-9   1e-3   3    1;  ...  % 10  wac
        ]; 
[fP,GAS,bi,bf]= gaoptimize(@MMC_Fitness_v3,GAP,D);

%% save the design
save('Design_5_ESTS_5kV', 'fP','GAS','bi','bf','D','GAP')

% [fP,GAS,bi,bf]= gaoptimize(@MMC_Fitness_v3,GAP,D, GASi, iP);
% save('Design_4_ESTS_5kV', 'fP','GAS','bi','bf','D','GAP')
