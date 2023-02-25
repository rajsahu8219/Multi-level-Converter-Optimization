function [f] = MMC_Fitness_v3(p,D,fn)
% This fundtion does MMC optimization, calling metamodel and EI core design
% code to optimize the MMC

% List of parameters
% First Set of parameters belong to the MMC

%  p(1) : Minimum Linc from MEC
%  p(2) : Submodule Capacitance
%  p(3) : Number of switching cycles per fundamental cycle
%  p(4) : Number of Submodules
%  p(5) : DC Bus Capacitance
%  p(6) : DC Bus Voltage -- variable for constant power

% Parameters for EI Core Inductor
% Will be lumped as P
%  p(7)    = Core material specification (an integer)
%  p(8)    = Conductor material specification (an integer)
%  p(9)    = Center leg width (m)
%  p(10)   = Twice end leg width to center leg width ratio
%  p(11)   = Twice I core width to center leg width ratio
%  p(12)   = Twice base of E core width to center leg width ratio
%  p(13)   = Desired conductor cross sectional area (m^2)
%  p(14)   = Desired number of turns
%  p(15)   = Desired number of conductors in width
%  p(16)   = Desired number of conductors in height
%  p(17)   = length of permanent magnet material (m)
%  p(18)   = Material type of permanent magnet
%  p(19)   = Length of the core (m)
%  p(20)   = air gap length (m)
%  p(21)   = dw/ds
%  p(22)   = ww/ws;

% Parameters for generator
% p(23) : Machine Speed in RPM
% p(24) : Desired Machine Specific Torque Density 

%% Open the parameters vector
% Generator
wrpm      = round(p(1)); % Speed of the machine in RPM
roTM      = p(2); % Desired Specific Torque Density of machine
Plmmc_est = p(3); % estimate MMC losses
Ns1a      = p(4); % 

% MMC
Csm      = p(5);
Nsm      = round(p(6));
Nscpfc   = round(p(7));
Lsi      = p(8);

% Bus capacitance
Cbus     = p(9);

% Inductor Metamodel
wac = p(10);

% MMC Modulation
ModType = 0;

% Constraints
CI = 0; CS = 0;
NC = 18;

% -------------------------------------------------------------------------

%% First, interact with Metamodel to find a machine for us

Pdc = D.Pdc; % this is the DC power required from us
Vdc = D.Vdc;
% Equivalent resistance
Rload = Vdc^2/Pdc;

% Generator power output
Pgo = Plmmc_est + Pdc;

% This optimization only considers electromagnetic mass, volume and loss of the
% machine
SMPMG = D.SMPMG;

% Feed DC Bus Voltage 
SMPMG.Vdc = Vdc;

% Variable Ns1a
SMPMG.Ns1a = Ns1a;

% convert speed to rad/s
wrmr = wrpm*(2*pi/60);

% call the electric machine metamodel 
[MPm,~,PrM,~,~,PlM,etaM,Jm,~,~,tipsp,Rsm,Lqm,Ldm,lambdam,~,PolePair,~] = SMPM_generator_metamodel_v1(SMPMG,roTM,Pgo,wrmr);

% check if any of the number is not a number and put a constant
if isnan(MPm) || MPm<0 || PrM<0 || PlM <0 || isnan(PlM) || isnan(etaM) || isnan(Rsm) ...
        || isnan(Lqm) || isnan(Ldm) || isnan(lambdam) || isnan(PolePair)
    cg1 = 0;
else
    cg1 = 1;
end

% frequency
PolePair = round(PolePair);

% Extract useful data about machine
% Electrical speed
wre = wrmr*PolePair;
fre = wre/(2*pi);

CS = CS+ cg1 + lte(fre, D.fundlim) + lte(Jm, D.Jpeak_m)+ lte(tipsp, D.tipspmax);
CI = CI+4;
if CS<CI
    % go out of the solving
     f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end

% Machine q-axis current
Iqsr = -(Pgo+PrM)/(1.5*wre*lambdam);

% Machine Current  -- assumed
Idsr     = 0;
Iqd = [Iqsr; Idsr];
Pin_est = Pgo+PlM;
% Pin1 = Pgo+PlM;

% % temp
% Iqtemp = sqrt(2)*sqrt((PrM/3)/(Rsm));
% tic

%% Next, find Inductor Metamodel Parameters Assuming major DC and fundamental component
% Calculate DC load Current
Idc = Pdc/Vdc;

% Define Machine Current
Is = sqrt(Iqsr^2 + Idsr^2);
phii1 = atan2(-Idsr, Iqsr); % phase angle of a phase
%% Inductor Calculations
% Following part is copy pasted from EI_fit_V1
% Do these calculations for only one of the phases and one of the legs
% otherwise following has to be done 6 times and we understand from power
% balance that the single phase inductors will observe near equal power
% loss
% I am picking up phase a upper leg -- any other phase leg could be picked
% -------------------------------------------------------------------------

% approximation mass calculation
EmL = (1/2)*Lsi*((Idc/3 +Is/2))^2;
MEI = D.L.cM*EmL*prod((D.L.bM + ( (Idc/3 +Is/2)/wac )*EmL^(1/3)).^(D.L.nM));


Kj = (  sqrt( (Idc/3)^2 +  ((Is/2))^2/2  )  )/(Idc/3 +Is/2);

PEIRLoss = D.L.cP*Kj^2*EmL^(1/3)...
    *prod((D.L.bP + ((Idc/3 +Is/2)/wac)*EmL^(1/3)).^(D.L.nP));

RLs = PEIRLoss/(( (Idc/3)^2 +  ((Is/2))^2/2  )  );
%% Next, do MMC simulation for an assumed L
D.L.R1 = RLs; % DC Resistance of the indcutor coil
W.Rcl  = D.L.R1;

% % Find equaivalent R load 
Rdc = Pdc/Vdc;

% Find DC voltage stress per SM--------------------------------------------
VperSm = Vdc/Nsm;

% Extract switch voltage rating of switch
Vlevel = D.sw_data(:,1);
temp1  = find(Vlevel >= VperSm);
if isempty(temp1)
    lvl = length(Vlevel);
else
    lvl = temp1(1);
end

% Conduction Loss Parameters
D.M.vfsw = D.sw_data(lvl, 7);
D.M.rsw  = D.sw_data(lvl, 8);
D.M.vfd  = D.sw_data(lvl, 9);
D.M.rd   = D.sw_data(lvl, 10);

% New Conduction Loss Parameters
% [beta; alpha; gamma]
% [Switch Diode]
tempconddata = D.swdata2(lvl, :);
% vdiode_up   = D.M.d.b + D.M.d.a.*iupapp .^D.M.d.n;
% vswitch_up  = D.M.s.b + D.M.s.a.*iupapp .^D.M.s.n;
% vdiode_low  = D.M.d.b + D.M.d.a.*ilowapp.^D.M.d.n;
% vswitch_low = D.M.s.b + D.M.s.a.*ilowapp.^D.M.s.n;
D.M.s.b = tempconddata(1); D.M.s.a = tempconddata(2); D.M.s.n = tempconddata(3);
D.M.d.b = tempconddata(4); D.M.d.a = tempconddata(5); D.M.d.n = tempconddata(6);

% Switching loss parameters
% normalized loss
Eofn = D.sw_data(lvl,4);
Eonn = D.sw_data(lvl,5);
Ern  = D.sw_data(lvl,6);

% temperature -- thermal resistances
% Rthsw2c = D.sw_data(lvl,11)+D.sw_data(lvl,13);
% Rthd2c  = D.sw_data(lvl,12)+D.sw_data(lvl,13);
Rthsw2c = (D.sw_data(lvl,11) + D.sw_data(lvl,13))/2; % quarter as 2 devices are coupled together
Rthd2c  = (D.sw_data(lvl,12) + D.sw_data(lvl,13))/2;

% Each Submodule gate driver mass
Msmpgd = (D.sw_data(lvl,14) + D.sw_data(lvl,15))*4; % 4 devices in a sm

% find switching frequcny --------------------------------------------
fsw = (wre/(2*pi))*Nscpfc;
% fprintf('Switching Freq : %d kHz \n',fsw/1e3)

% establish number of points to represent per cycle,--------------------------------------------
% time, angle of synchronous reference frame, fft matrix
Np = Nscpfc*D.N.Nppsw;
t  = linspace(0,2*pi/wre,Np+1);
t  = t(1:end-1);

% find angles, Transformation matrices----------------------------------
qe0 = wre*t;
qeM = qe0 - 2*pi/3;
qeP = qe0 + 2*pi/3;

theta_vec  = [qe0;qeM;qeP];

% Find Current Variables --------------------------------------------
iabcs = Is*cos(theta_vec + phii1);

% estimate the upper leg current
iupapp  = -Idc/3 + iabcs/2 ;

% estimate the lower leg current
ilowapp = -Idc/3 - iabcs/2 ;

% Find the remaining parameters for Calculations--------------------------------------------
% Inductance Paramenters
Lmi = 0; % zero mutual inductance
Lleg = [Lsi Lmi Lmi; Lmi Lsi Lmi; Lmi Lmi Lsi];
Lqeq = Lqm + (Lsi-Lmi)/2;
Ldeq = Ldm + (Lsi-Lmi)/2;

% Resistance Parameters - overestimated to calculate the duty cycle
% equaivalent switch resistance
Rsw  = 2*Nsm*max(D.M.rsw,D.M.rd); 
Vmd  = 2*Nsm*max(D.M.vfsw,D.M.vfd); 
Rmmc = Rsw + W.Rcl;
Req  = Rsm + Rmmc/2; % removed resistance in inductor

% Find the Machine Side Small Signal Impedence--------------------------------------------
Zeq  = [Req, wre*Ldeq;  -wre*Lqeq, Req];
% Ziqdac = bsxfun(@plus, bsxfun(@times, s1, [Lqeq , 0; 0, Ldeq]), Zeq );


% find refernece frame : assuming synchronous refernce frame
Ksinvfc    = [cos(theta_vec)];
Ksinvsc    = [sin(theta_vec)];

% compute the triangle wave--------------------------------------------
qsw  = qe0*Nscpfc;
wup  = acos(cos(qsw))/pi;
% wlow = acos(cos(qsw+pi))/pi;
wlow = wup; % temp

% Find shift in Triangles
trishift = linspace(Nsm-1,0,Nsm)';

% find total triangle
triup  = bsxfun(@plus, wup, trishift)/Nsm;
trilow = bsxfun(@plus, wlow, trishift)/Nsm;

% MOdulation Type -- Store it in D 



% find the angles where zero crossing happeS.Ns
% first check if the machine current makes sense
% constrain 2: -check if idc/(3*iac1) is between -1 and 1
cm1 = lte(abs(Idc/(3*Is/2)),1);
CS = CS+cm1;
CI = CI+1;
if CS<CI
    % go out of the solving
    f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end

zc = acos(Idc/(3*Is/2));

% Find Fundamental Component of Converter Voltage -------------------------
Vconvqdr = Zeq*Iqd + [wre*lambdam; 0];

% Approximate Capacitor Voltages --------------------------------------------
% Find the harmonic components in form of ac1*cos(theta) - bc1*sin(theta)
% estimate the capacitor voltages
% voltage ripple fundamental component
abc1 = -[Iqd(2); Iqd(1)]/(4*wre*Csm) - [Vconvqdr(2); Vconvqdr(1)]*Idc/(3*wre*Csm*Vdc);

% voltage ripple Second harmnonic component
abc2 = [Iqd(2), Iqd(1); Iqd(1), -Iqd(2)]*Vconvqdr *1/(8*Vdc*wre*Csm) ;

% find DC component of leg voltage -- Note to self capacitor current does
% not have a DC component
% assume the voltage of the voltage kVdc/N
k1 = 1 + Nsm/(Vdc)^2*[abc1(1), -abc1(2)]*Vconvqdr - 2*Vmd*(2*zc/pi-1)/Vdc +  2*(Rmmc)/(3*Rload);
% k1 = 1 + Nsm/(Vdc)^2*[abc1(1), -abc1(2)]*Vconvqdr ;

% Find second harmonic compoenent -----------------------------------------
vr2a2 = Nsm/2*(abc2(1) - abc1(1)*Vconvqdr(1)/Vdc - abc1(2)*Vconvqdr(2)/Vdc)/k1 + Vmd*cos(2*phii1)*( 2*sin(2*zc)/pi)/k1;
vr2b2 = Nsm/2*(abc2(2) + abc1(1)*Vconvqdr(2)/Vdc - abc1(2)*Vconvqdr(1)/Vdc)/k1 + Vmd*sin(2*phii1)*( 2*sin(2*zc)/pi)/k1;

% a phase converter voltage
vconv    = Vconvqdr(1)*Ksinvfc + Vconvqdr(2)*Ksinvsc;

% Estimate the second harmonic component to be removed
v2conv = vr2a2*cos(2*theta_vec) - vr2b2*sin(2*theta_vec);

% Esimate DC + Fundamental Refernece Voltage ------------------------------
% find a phase upper leg reference submodule voltage and duty cycle
% vupref  = bsxfun(@plus, -vconv - v2conv, Vdc/2);
vupref  = bsxfun(@plus, -vconv - v2conv, Vdc/2);

% find a phase lower leg reference submodule voltage
% vlowref = bsxfun(@plus, vconv - v2conv, Vdc/2);
vlowref = bsxfun(@plus,  vconv - v2conv, Vdc/2);

% find modulation index
mup = vupref/(Vdc);
mlow = vlowref/(Vdc);

% limit the maximum of modulation index less than 1
m_max = max(max(abs(mup(:))),max(abs(mlow(:))));

% constraint 2:  check if mnodulation index is <=1
cm2 = lte(m_max,1);
% constraint 3 : check if fsw is less than specfied
cm3 = lte(fsw, D.fswmax);
CS = CS + cm2+cm3;
CI = CI+2;
if CS<CI
    % go out of the solving
    f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end

% Calculate Sum of Switching Signals
% Calculate total voltage across submodules
% for upper leg switches
% Calculate total voltage across submodules
% for upper leg switches
swup = zeros(size(mup)); 
swlow = swup;
% swup = cell(3,Nsm); swlow = cell(3,Nsm);
% sm_aswup = 0;
for i_it = 1:Nsm
    swup(1,:) =  swup(1,:) + ones(size(mup(1,:))).*(mup(1,:)>=triup(i_it,:));
    swup(2,:) =  swup(2,:) + ones(size(mup(1,:))).*(mup(2,:)>=triup(i_it,:));
    swup(3,:) =  swup(3,:) + ones(size(mup(1,:))).*(mup(3,:)>=triup(i_it,:));
%     sm_aswup = swup{1,i_it}+sm_aswup; % bug fixes - needed to calculate cap I
end

% Find Switch Signals for lower leg 

if ModType ==1
    
    for i_it = 1:Nsm
        swlow(1,:) = swlow(1,:) + ones(size(mup(1,:))).*(mlow(1,:)>=trilow(i_it,:));
        swlow(2,:) = swlow(2,:) + ones(size(mup(1,:))).*(mlow(2,:)>=trilow(i_it,:));
        swlow(3,:) = swlow(3,:) + ones(size(mup(1,:))).*(mlow(3,:)>=trilow(i_it,:));
    end
    
else
    
    for i_it = 1:Nsm
        swlow(1,:) = swlow(1,:) + ones(size(mup(1,:))).*(mlow(1,:)>=triup(i_it,:));
        swlow(2,:) = swlow(2,:) + ones(size(mup(1,:))).*(mlow(2,:)>=triup(i_it,:));
        swlow(3,:) = swlow(3,:) + ones(size(mup(1,:))).*(mlow(3,:)>=triup(i_it,:));
    end
    
end

% estimate the upper leg capacitor voltage
vcapup  = k1*Vdc/Nsm + abc1(1)*cos(theta_vec) - abc1(2)*sin(theta_vec) + abc2(1)*cos(2*theta_vec) - abc2(2)*sin(2*theta_vec) ;

% estimate the lower leg voltage
vcaplow = k1*Vdc/Nsm - abc1(1)*cos(theta_vec) + abc1(2)*sin(theta_vec) + abc2(1)*cos(2*theta_vec) - abc2(2)*sin(2*theta_vec) ;

% estimating capacitor voltage ripples in upper and lower leg --NOTICE
vriplow  = max(vcaplow(1,:)) - min(vcaplow(1,:));
vcmaxlow = max(abs(vcaplow(1,:)));

% find capcitor voltage ripple and max cap voltage for phase a -- assume it same for all the phases
vripup  = max(vcapup(1,:)) - min(vcapup(1,:));
vcmaxup = max(abs(vcapup(1,:)));

vcripmax = max(vripup,vriplow);
vcmax    = max(vcmaxup,vcmaxlow);

% constraint 4:  limit max capacitor voltage ripple
cm4 = lte(vcripmax,(D.C.vrip)*Vdc/Nsm);
CS  = CS + cm4;
CI  = CI+1;
if CS<CI
    % go out of the solving
    f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end

% MMC Leg Voltages  --------------------------------------------

% Semiconductor device drops
vdiode_up   = D.M.d.b + D.M.d.a.*abs(iupapp) .^D.M.d.n;
vswitch_up  = D.M.s.b + D.M.s.a.*abs(iupapp) .^D.M.s.n;
vdiode_low  = D.M.d.b + D.M.d.a.*abs(ilowapp).^D.M.d.n;
vswitch_low = D.M.s.b + D.M.s.a.*abs(ilowapp).^D.M.s.n;

% Sqwaveup = Sqwave*((iabcs/2 - D.Idc/3)>0)*-1;
NullMat =  zeros(size(iabcs));
vswup_diode = (iupapp >= NullMat).*(swup.*2.*vdiode_up + (Nsm*ones(size(swup)) -swup ).*(vdiode_up ) )...
    + (-1)*(iupapp < NullMat).*( (Nsm*ones(size(swup)) -swup ).*(vdiode_up ) ) ;
vswup_switch  = (iupapp >= NullMat).*( (Nsm*ones(size(swup)) -swup ).*( vswitch_up) )...
    + (-1)*(iupapp < NullMat).*(swup.*2.*vswitch_up + (Nsm*ones(size(swup)) -swup ).*( vswitch_up) ) ;
% vswup  = (iupapp >= NullMat).*(swup.*2.*vdiode_up + (Nsm*ones(size(swup)) -swup ).*(vdiode_up + vswitch_up) )...
%     + (-1)*(iupapp < NullMat).*(swup.*2.*vswitch_up + (Nsm*ones(size(swup)) -swup ).*(vdiode_up + vswitch_up) ) ;

% separate for heat sink calculatins
%-- Diode upper leg
% find diode voltage drop when SM was ON and i>0
vswup_diode_iplus_smon = (iupapp >= NullMat).*(swup.*2.*vdiode_up);
% find diode voltage drop when SM was OFF and i>0
vswup_diode_iplus_smoff = (iupapp >= NullMat).*((Nsm*ones(size(swup)) -swup ).*(vdiode_up ));
% find diode voltage drop when SM was OFF and i<0
vswup_diode_iminus_smoff = (-1)*(iupapp < NullMat).*( (Nsm*ones(size(swup)) -swup ).*(vdiode_up ) );
% fourth case does not exist for diode

%-- switch upper leg
% find switch voltage drop when SM was ON and i<0
vswup_sw_iminus_smon = (-1)*(iupapp < NullMat).*(swup.*2.*vswitch_up);
% find diode voltage drop when SM was OFF and i<0
vswup_sw_iminus_smoff = (-1)*(iupapp < NullMat).*( (Nsm*ones(size(swup)) -swup ).*( vswitch_up) );
% find switch voltage drop when SM was OFF and i>0
vswup_sw_iplus_smoff = (iupapp >= NullMat).*( (Nsm*ones(size(swup)) -swup ).*( vswitch_up) );
% fourth case does not exist for diode

vswlow_diode = (ilowapp >= NullMat).*(swlow.*2.*vdiode_low + (Nsm*ones(size(swlow)) -swlow ).*(vdiode_low ) )...
    +(-1)* (ilowapp < NullMat).*((Nsm*ones(size(swlow)) -swlow ).*(vdiode_low ) ) ;

vswlow_switch = (ilowapp >= NullMat).*( (Nsm*ones(size(swlow)) -swlow ).*( vswitch_low) )...
    +(-1)* (ilowapp < NullMat).*(swlow.*2.*vswitch_low + (Nsm*ones(size(swlow)) -swlow ).*( vswitch_low) ) ;

% vswlow = (ilowapp >= NullMat).*(swlow.*2.*vdiode_low + (Nsm*ones(size(swlow)) -swlow ).*(vdiode_low + vswitch_low) )...
%     +(-1)* (ilowapp < NullMat).*(swlow.*2.*vswitch_low + (Nsm*ones(size(swlow)) -swlow ).*(vdiode_low + vswitch_low) ) ;
% vswlow  = Vmd*(2*(( ilowapp ) >= NullMat) - 1);

% Calculate total voltage across submodules
% for upper leg switches

vup = bsxfun(@times,   vcapup, [swup(1,:); swup(2,:); swup(3,:)]);


% finding equivalent voltage drop in devices

% % vaup = vcapup.*Nsm.*maup;
vup  = vup + vswup_diode + vswup_switch + W.Rcl*iupapp;

% for lower leg switches
vlow =  bsxfun(@times, vcaplow, [swlow(1,:); swlow(2,:); swlow(3,:)]);


% % valow = vcaplow.*Nsm.*malow;

vlow = vlow + vswlow_diode + vswlow_switch + W.Rcl*ilowapp;

% perform the reduced order simulation to estimate the switching freq leg currents
U = bsxfun(@minus, Vdc, vup+vlow);
U = U - mean(U,2); % get rid of any numerical issues

% Check if the flag is true
% if flag1==1 % most likely the zero sequence impedance is zero
%     temp = NaN*size(U);
% else
dtemp = Lleg\U;
temp  = cumtrapz(t, dtemp,2);
DC_bias = mean(temp,2);
temp =  temp - DC_bias+ (-2*Idc/3);
% end

iabcup  = (iabcs+temp)/2;
iabclow = (temp-iabcs)/2;

% Impose constraint on maximum current ripple? 
var1 = max(abs(iabclow(1,:)-mean(iabclow(1,:)))); var2 = max(abs(iabcup(1,:)-mean(iabcup(1,:))));
cm5   = lte(max(var1,var2),D.irip*Is/2);

% Put constraint on power imbalanbce
cm6 = lte(max(abs(DC_bias)), 0.01*abs(2*Idc/3));

CS = CS + cm5 + cm6;
CI = CI+2;

% if CS<CI
%     % go out of the solving
%     f = 1e-2*[1 1 1]*(CS-NC)/NC;
%     return
% end

% calculate MMC Conduction losses  --------------------------------------------
PlossMMC3phase = fre*trapz(t,(vswup_diode + vswup_switch).*iupapp, 2 ) +  fre*trapz(t,(vswlow_diode + vswlow_switch).*ilowapp,2 ) ;
Plmmc_cond = sum(PlossMMC3phase);

% put another constraint on power imbalaance
Plossdiff = max(abs(diff(PlossMMC3phase)))/mean(PlossMMC3phase);


% One last constraint on DC bus capacitance
idcbus = -sum(iabcup,1)-Idc;

vdcbus_est = (cumtrapz(t,idcbus))/Cbus + Vdc;
vdcrip_est = max(vdcbus_est) - min(vdcbus_est);

cm7  =  lte(vdcrip_est, D.vbusrip*Vdc);

% Impose constraint on Maximum fault ridethrough current %%%% future

CS = CS + cm7;
CI = CI+1;

if CS<CI
    % go out of the solving
    f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end
% toc
%% Semiconductor switching losses
% Switching loss -  approximation in phase a upper leg
% Calculation for upper leg
Eloss_1leg_sw = squeeze((Eonn + Eofn)*(Vdc/Nsm)*fre*( trapz(t,abs(iabcup(1,:))) )); 
Eloss_1leg_d = squeeze(( Ern)*(Vdc/Nsm)*fre*( trapz(t,abs(iabcup(1,:))) ));
PSw_loss = fsw*(Eloss_1leg_sw+Eloss_1leg_d);
%% Polyprop Capacitor Calculations
% Mass Calculation  --------------------------------------------
Mcap = D.C.M.alpha*Csm.*(((Vdc/Nsm)./D.C.M.vb).^D.C.M.gamma) + D.C.M.beta;
% Mcap = D.C.betac*Csm*(vcmax)^1.5;

% estimate actual capacitor current ripple: consider only the fundamental
% component of the AC side current
% for one of the capacitors
iacapup    = swup(1,:).*squeeze(iabcup(1,:))/Nsm;
iacapuprms =  sqrt(fre*(trapz(t, iacapup.^2 )));%rms(iacapup);

% Losses in capacitor - Neglecting its effect on MMC branch drop  ---------
% using the model to calculate capacitor resistance D.C.ESR.
Rcap = (D.C.ESR.alpha1*(D.C.ESR.vb1./(Vdc/Nsm)).^D.C.ESR.gamma1 + D.C.ESR.alpha2*(D.C.ESR.vb2./(Vdc/Nsm)).^D.C.ESR.gamma2)./Csm + D.C.ESR.beta;
% Rcap = D.C.R.alpha/((Csm^D.C.R.n1)*(Vdc/Nsm)^D.C.R.n2); % corrected
Pcaprloss = Rcap*Nsm*iacapuprms^2;

% heat conduction in capacitors
% Thermal Conductivity Calculation
G_csm = D.C.agc*((Vdc/Nsm)/D.C.Gvbase).*(Csm).^D.C.aga + D.C.agd*((Vdc/Nsm)/D.C.Gvbase).*(Csm).^D.C.agb ...
    +   D.C.bgc*(D.C.Gvbase/(Vdc/Nsm)).*(Csm).^D.C.bga + D.C.bgd*(D.C.Gvbase/(Vdc/Nsm)).*(Csm).^D.C.bgb; 
Delta_T_Csm = (Pcaprloss/(Nsm*G_csm));

% DC bus losses ND mass
Nbcap = ceil(Vdc/D.C.Vcdsmax); % number of DC bus capacitors

Rcdc = Nbcap*(( D.C.ESR.alpha1*(D.C.ESR.vb1./(Vdc/Nbcap)).^D.C.ESR.gamma1 + D.C.ESR.alpha2*(D.C.ESR.vb2./(Vdc/Nbcap)).^D.C.ESR.gamma2 )./(Nbcap*Cbus) + D.C.ESR.beta );
Pcdcbusloss = Rcdc*(fre*(trapz(t, idcbus.^2 )));  % fix resistive losses
Mcdc = Nbcap*(D.C.M.alpha*(Nbcap*Cbus).*(((Vdc/Nbcap)./D.C.M.vb).^D.C.M.gamma) + D.C.M.beta);

% heat conduction in capacitors

G_cin = D.C.agc*((Vdc/Nbcap)/D.C.Gvbase).*(Cbus).^D.C.aga +D.C.agd*((Vdc/Nbcap)/D.C.Gvbase).*(Cbus).^D.C.aga ...
     +  D.C.bgc*(D.C.Gvbase/(Vdc/Nbcap)).*(Cbus).^D.C.bga +D.C.bgd*(D.C.Gvbase/(Vdc/Nbcap)).*(Cbus).^D.C.bga; 

Delta_T_Cin = (Pcdcbusloss/(Nbcap*G_cin));

% constraints
cc1 = lte(Delta_T_Csm, D.C.DTempmax);
cc2 = lte(Delta_T_Cin, D.C.DTempmax);
%% Remaining Inductor Calculations

ilowmax = max(abs(iabclow(1,:))); iupmax = max(abs(iabcup(1,:))); 
ilegmax = max(ilowmax,iupmax);
Jpk = (ilegmax/wac);

ce1 = lte(Jpk,D.Jpeak); 



CS = CS + cc1+ cc2+ ce1;
CI = CI+3;

if CS<CI
    % go out of the solving
    f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end


% disp('Went Through Indcutor')
%--------------------------------------------------------------------------
%% Heat Sink --------------------------------------------------------------
% Average Conduction loss calculations -- Diode pairs 
% D1 plus D4
PcLd1d4 = (sum(mean(vswup_diode_iplus_smon.*iupapp,2))+sum(mean(vswup_diode_iplus_smoff.*iupapp,2)))/(3*Nsm);
% D2 plus D3
PcLd2d3 = (sum(mean(vswup_diode_iminus_smoff.*iupapp,2)))/(3*Nsm);
% Average Conduction loss calculations -- IGBT pairs 
% S1 plus S4
PcLs1s4 = (sum(mean(vswup_sw_iminus_smon.*iupapp,2))+sum(mean(vswup_sw_iminus_smoff.*iupapp,2)))/(3*Nsm);
% S2 plus S3
PcLs2s3 = (sum(mean(vswup_sw_iplus_smoff.*iupapp,2)))/(3*Nsm);

% Average switching loss
% D1 plus D4
Psloss_d_iplus   = (fsw/Nsm)*squeeze( Ern*(Vdc/Nsm)*mean(abs(iabcup(1,:).*((iabcup(1,:)>=0)))));
% D2 plus D3
Psloss_d_iminus  = (fsw/Nsm)*squeeze( Ern*(Vdc/Nsm)*mean(abs(iabcup(1,:).*((iabcup(1,:)<0)))));
%S1 plus S4
Psloss_sw_iminus = (fsw/Nsm)*squeeze((Eonn + Eofn)*(Vdc/Nsm)*mean(abs(iabcup(1,:).*((iabcup(1,:)<0)))));
%S2 plus S3
Psloss_sw_iplus  = (fsw/Nsm)*squeeze((Eonn + Eofn)*(Vdc/Nsm)*mean(abs(iabcup(1,:).*((iabcup(1,:)>=0)))));

% Total average losses in device pairs in a submodule
Pd1d4 = PcLd1d4+Psloss_d_iplus;
Pd2d3 = PcLd2d3+Psloss_d_iminus;
Ps1s4 = PcLs1s4+Psloss_sw_iminus;
Ps2s3 = PcLs2s3+Psloss_sw_iplus;

% Rthsw2c Rthd2c
% Figure out heat sink for junction temperature of device pairs
Rth = (1./(Pd1d4+Pd2d3+Ps1s4+Ps2s3))*( D.HS.Tjmax - D.HS.Tamb - [Pd1d4*Rthd2c; Pd2d3*Rthd2c; Ps1s4*Rthsw2c; Ps2s3*Rthsw2c] );

MinRth = min(Rth);

if MinRth>D.HS.Rlim % if not a lot of heat sink is required, use the least mass possible
    MinRth = D.HS.Rlim;
end
chs1 = gte(MinRth, D.HS.Rmin); % to make sure heat transfer is feasible

CS = CS + chs1   ;
CI = CI + 1;
if CS<CI
    % go out of the solving
   f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end

MaxMhs =  (D.HS.a./(MinRth).^D.HS.n1 + D.HS.b./(MinRth).^D.HS.n2);
Pfan = Nsm*(D.HS.f.a./(MinRth).^D.HS.f.n1 + D.HS.f.b./(MinRth).^D.HS.f.n2);

chs2 = lte(MaxMhs, D.HS.Mmax);

%% total Switch and gate drive mass
Msw = Msmpgd*6*Nsm;
%% Total Loss and Mass-----------------------------------------------------
% Does include switching losses--------------------------------------------
Loss_total =  PlM + (Plmmc_cond + 6*PSw_loss + 6*Pfan + 6*PEIRLoss + 6*Pcaprloss) + Pcdcbusloss;
Plmmc_act = (Plmmc_cond + 6*PSw_loss + 6*Pfan + 6*PEIRLoss + 6*Pcaprloss) + Pcdcbusloss;

% Input Powers
Pin_act = Loss_total + Pdc; % (actual)

% Efficiency
Eta = 100*(Pdc)/(Pin_act);

% Total Mass -----------------------------------------------------
Mass_MMC = 6*MEI*D.MfacI + 6*Nsm*(Mcap+ MaxMhs+Msmpgd)*D.MfacSM + Mcdc;
Mass_total = MPm + Mass_MMC;


% Constraint on MMC solution convergence
% Constraint for a viable design which satisfies power requirements
cp1 = gte(Pin_est,Pin_act);

cp2 = lte(Mass_total, D.Mmax);


CS = CS + chs2 + cp1 + cp2 ;
CI = CI+3;

if (NC < CI) || (NC>CI)
    disp(['NC = ' num2str(NC)])
    disp(['CI = ' num2str(CI)])
    error('Mismatch in constraints')
end
if CS<CI
    % go out of the solving
    f = 1e-2*D.UnitMetric*(CS-NC)/NC;
    return
end

% Power input avergae metric
Loss_metric = Loss_total + 0.5*( Plmmc_est - Plmmc_act )^2;

%% objective function formulation
if (numel(D.UnitMetric))>2
    f = [1/(Mass_total) 1/Loss_metric  1/(MEE) 1/(P_EE) ];
else
    f = [1/(Mass_total) 1/Loss_metric ];
end


%% Post processing
if nargin>2
    
    clear f
    
    % Mass
    f.Mt    = Mass_total;
    f.Mpm   = MPm;
    f.MEI   = MEI*D.MfacI;
    f.Mcap  = Mcap;
    f.Mcdc = Mcdc;
    f.Mhs   = MaxMhs;
    f.Mswp  = Msmpgd; % mass of one switch package
    f.Mmmc  = Mass_MMC;
    f.MSM   = (Mcap+ MaxMhs+Msmpgd)*D.MfacSM;
    
    % Loss
    f.Eta   = Eta;
    f.Lt    = Loss_total;
    f.Lmc   = PlM;
    f.LMMC  = Plmmc_cond;
    f.LSw   = PSw_loss;
    f.LEIr  = PEIRLoss;
%     f.LEIr  = PEIrloss;
    f.LCap  = Pcaprloss;
    f.LDCbCap = Pcdcbusloss;
    
    % Machine Parameters
    f.PM.Rsm = Rsm;
    f.PM.Lqm = Lqm;
    f.PM.Ldm = Ldm;
    f.PM.lambdam = lambdam;
    f.PM.PolePair = PolePair;
    f.PM.wrmrpm = round(wrpm);
    f.Iqsr = Iqsr;
    f.Poutac = Pgo;
    f.Pin = Pin_est;
    f.Pin_act = Pin_act;
    f.Jm = Jm;
    
    % Inductor Parameters
%     f.EI   = EI;
%     f.W    = W;
%     f.Linc = Linc;
    f.wac  = wac;
    
    % MMC Parameters
    f.Pdc = Pdc;
    f.Csm = Csm;
    f.Cdc = Cbus;
    f.Nsm = Nsm;
    f.Lsind = Lsi;
    f.fsw = fsw;
    f.Vdc = Vdc;
    
    
    if fn>0
        
        % plot MMC plots
        
        figure(fn)
        plot(t,iabcs)
        title('Machine Current vs Time')
        xlabel('Time (s)')
        ylabel('Current (A)')
        fn = fn+1;
        
        figure(fn)
        plot(t,iabcup)
        title('MMC Upper Leg Current vs Time')
        xlabel('Time (s)')
        ylabel('Current (A)')
        fn = fn+1;
        
        figure(fn)
        plot(t,iabclow)
        title('MMC Lower Leg Current vs Time')
        xlabel('Time (s)')
        ylabel('Current (A)')
        fn = fn+1;
        
        figure(fn)
        plot(t,vcapup)
        title('Upper Leg Avg. Capacitor vs Time')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        fn = fn+1;
        
        figure(fn)
        plot(t,vcaplow)
        title('Lower Leg Avg. Capacitor vs Time')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        fn = fn+1;
        save('Data_52_1MW','iabcup','iabcs', 'vcapup','t')
        
%         figure(fn)
%         plot(Linc/Lsind)
%         title('Ratio of Increamental and Assumed Inductance vs Time')
%         xlabel('Time (s)')
%         ylabel('Ratio of Increamental and Assumed Inductance')
%         fn = fn+1;
%         
%         draw_EI_PM(EI,W,2,fn)
        
        
        disp([   ]);
        disp('---------------------------- Design Data ----------------------------');
        disp([' System Mass ' num2str(Mass_total) ' kg '])
        disp([' System Loss ' num2str(Loss_total/1e3) ' kW '])
        disp('---------------------------- Generator Data ----------------------------');
        disp([' Generator Speed ' num2str(round(wrpm)) ' RPM '])
%         disp([' Generator Power Rating ' num2str(Poutest/1e3) ' kW '])
        disp([' Generator Torque Density ' num2str(roTM) ' Nm/kg ' ])
        disp(['      Generator Mass ' num2str(MPm) ' kg ']);
        disp(['      Generator Total Loss ' num2str(PlM/1e3) ' kW ']);
        disp(['      Generator Rs ' num2str(Rsm*1e3) ' mOhm ']);
        disp(['      Generator Lq ' num2str(Lqm*1e3) ' mH ']);
        disp(['      Generator lambdam ' num2str(lambdam*1e3) ' mVs ']);
        disp(['      Generator Pole Pairs ' num2str(PolePair) '  ']);
        disp([   ]);
        disp('---------------------------- MMC Data ----------------------------');
        disp([' MMC Estimated Leg inductance ' num2str(Lsi/1e-3) ' mH ']);
        disp([' Submodule capacitance ' num2str(Csm/1e-3) ' mF ']);
        disp([' DC Bus capacitance ' num2str(Cbus/1e-3) ' mF ']);
        disp([' Number of submodule ' num2str(Nsm)]);
        disp([' Switching frequency ' num2str(fsw/1e3) ' kHz ']);
        disp([' DC Bus Voltage ' num2str(Vdc/1e3) ' kV ']);
        disp(['      Total MMC Mass ' num2str(Mass_MMC) ' kg ']);
        disp(['      Individual Inductor Mass ' num2str( MEI*D.MfacI ) ' kg ']);
        disp(['      Individual Capacitor Mass ' num2str( Mcap ) ' kg ']);
        disp(['      Individual Heat Sink Mass ' num2str( MaxMhs ) ' kg ']);
        disp(['      Total MMC Loss ' num2str((Plmmc_cond + 6*PSw_loss + 6*PEIRLoss + 6*Pcaprloss)/1e3) ' kW ']);
        disp(['      MMC Leg Switch Conduction and Resitive Losses ' num2str((Plmmc_cond/6)/1e3) ' kW ']);
        disp(['      Each Leg Inductor Resistive Loss ' num2str(PEIRLoss/1e3) ' kW' ])
        disp(['      Leg Switching Loss ' num2str(( PSw_loss )/1e3) ' kW ']);
        disp(['      Leg Capacitor Loss ' num2str(( Pcaprloss )/1e3) ' kW ']);
      %  disp(['      Leg Indcutor Core Loss ' num2str(( PEIcloss )/1e3) ' kW ']);
        disp([   ]);
        disp('---------------------------- Inductor Data ----------------------------');
        
        disp(['    Wire cross sectional area: ' num2str(wac*10^6) ' mm^2']);
      %{  
        disp(['      Material = ' EI.MP.desc]);
        disp(['      wc (m) = ' num2str(EI.wc)]);
        disp(['      we (m) = ' num2str(EI.we)]);
        disp(['      wi (m) = ' num2str(EI.wi)]);
        disp(['      wb (m) = ' num2str(EI.wb)]);
        disp(['      ws (m) = ' num2str(EI.ws)]);
        disp(['      ds (m) = ' num2str(EI.ds)]);
        disp(['      lc (m) = ' num2str(EI.lc)]);
        disp(['      g (m) = ' num2str(EI.g)]);
        disp(['                           Winding Data']);
        disp(['      Material = ' W.MP.desc]);
        disp(['      AWG = ' W.desc]);
        disp(['      ac (m^2) = ' num2str(W.ac)]);
        disp(['      rc (m) = ' num2str(W.rc)]);
        disp(['      N = ' num2str(W.N)]);
        disp(['      Nw = ' num2str(W.Nw)]);
        disp(['      Nd = ' num2str(W.Nd)]);
        disp(['      dw (m) = ' num2str(W.dw)]);
        disp(['      ww (m) = ' num2str(W.ww)]);
        disp(['                          Metrics']);
     
        disp(['      Rcl (Ohms) = ' num2str(W.Rcl)]);
        disp(['      hE (m) = ' num2str(hE)]);
        disp(['      Lincmin = ' num2str(Lincmin)]);
        disp(['      wE (m) = ' num2str(wE)]);
        disp(['      lE (m) = ' num2str(lE)]);
        disp(['                    Permanent Magnet Data']);
        disp(['      Material = ' EI.PP.desc]);
        disp(['      Width of PM(m) = ' num2str(EI.lpm)]);
        disp([   ]);
        %}
         disp('---------------------------- Capacitor Data ----------------------------');
%         
%         disp([' Capacitor Life Time ' num2str(L2) ' Hours']);
%         disp([' Capacitor Current Ripple ' num2str(iacapuprms) ' A']);
        disp([' Capacitor Resistance ' num2str(Rcap) ' Ohm']);
        disp([' Capacitor Voltage Ripple ' num2str(vcripmax) ' V']);
        disp([   ]);
    end
    
end
end

function [c]=gte(x,xmnc)
%   FUNCTION:    [c]=gte(x,xmn)
%   DESCRIPTION: Function to test constraint. Returns 1 if the argument
%                is greater-than allowed value, otherwise it returns
%                a value between 0 and 1.
%   INPUTS:      x       -   variable to be measured
%                xmn     -   minimum allowed value for the variable
%   OUTPUTS:     f       -   fitness value

    if (x<xmnc)
        c=1/(xmnc-x+1);
    else
        c=1;
    end
     
end

function [c]=lte(x,xmxc)
%   FUNCTION:    [c]=lte(x,xmx)
%   DESCRIPTION: Function to test constraint. Returns 1 if the argument
%                is greater-than allowed value, otherwise it returns
%                a value between 0 and 1.
%   INPUTS:      x       -   variable to be measured
%                xmx     -   maximum allowed value for the variable
%   OUTPUTS:     f       -   fitness value

    if (x>xmxc)
        c=1/(x-xmxc+1);
    else
        c=1;
    end
end
