%DiMuon2017.m
%
% This program is being created to analyze CMS dimuon events
% GOAL: Find overall mass of Jpsi particle

clear; clc; close all;


%% READ THE EVENTS

%
% Read event number, muon quality and 4-momenta data from
% the Comma Separated Value (.csv) file
%

istart=5;                           % First row  is 0; therefore, sixth row is 5

Nevents=2000;                       % Number of events to read

EvtNum = csvread('IPPW Trinity 2017.csv',istart,1,[istart,1,istart+Nevents-1,1]);

MuQuality = csvread('IPPW Trinity 2017.csv',istart,2,[istart,2,istart+Nevents-1,2]);

P4Mu1 = csvread('IPPW Trinity 2017.csv',istart,3,[istart,3,istart+Nevents-1,6]);

P4Mu2 = csvread('IPPW Trinity 2017.csv',istart,7,[istart,7,istart+Nevents-1,10]);

%% PRINT THE 5TH EVENT TO CHECK THE READING OF THE DATA

% disp(EvtNum(5))
% disp(MuQuality(5))
% disp(P4Mu1(5,1:4))
% disp(P4Mu2(5,1:4))

%% INITIALIZE SOME ARRAYS TO ZERO

P4Mu1RF=zeros(2000,4);
P4Mu2RF=zeros(2000,4);
%% DIAGNOSTICS

%
% Diagnostic Plots: Histogram all 10 input variables
% in one Figure using subplots
%
figure('Name','Diagnostic Plot: 10 original input vars')
subplot(3,4,1)
histogram(EvtNum)
title('Event Numbers')

% subplot(3,4,2)

subplot(3,4,3)
histogram(MuQuality,4)
title('Quality of Muon Data')


% subplot(3,4,4)

subplot(3,4,5)
histogram(P4Mu1(:,1))
title('E1 Energy (GeV)')


subplot(3,4,6)
histogram(P4Mu1(:,2))
title('X-component for P1 (GeV/c)')

subplot(3,4,7)
histogram(P4Mu1(:,3))
title('Y-component for P1 (GeV/c)')

subplot(3,4,8)
histogram(P4Mu1(:,4))
title('Z-component for P1 (GeV/c)')

subplot(3,4,9)
histogram(P4Mu1(:,1))
title('E2 Energy (GeV)')

subplot(3,4,10)
histogram(P4Mu1(:,2))
title('X-component for P2 (GeV/c)')

subplot(3,4,11)
histogram(P4Mu1(:,3))
title('Y-component for P2 (GeV/c)')

subplot(3,4,12)
histogram(P4Mu2(:,4))
title('Z-component for P2 (GeV/c)')



%% DiMuon INVARIANT MASS - LOOP OVER THE EVENTS

    
    %
    % For each event calculate the DiMuon System's Energy (Esys),
    % 3-momentum (P3sys=[Px,Py,Pz]) and invariant mass
    %
Esys=P4Mu1(:,1)+P4Mu2(:,1);
P3sys=[P4Mu1(:,2)+P4Mu2(:,2)  ,  P4Mu1(:,3)+P4Mu2(:,3) ,  P4Mu1(:,4)+P4Mu2(:,4)];
Minv=sqrt(Esys.^2-sum(P3sys.^2,2));
Psys=sqrt(P3sys(:,1).^2+P3sys(:,2).^2+P3sys(:,3).^2);

% Histogram Esys, Psys (the magnitude of the 3-momentum)and the DiMuon
% invariant mass in one Figure using subplots
%
figure
subplot(3,1,1)
histogram(Esys)
title('Energy of the System (GeV)')
subplot(3,1,2)
histogram(Psys)
title('Magnitude of the 3-momentum (GeV/c)')
subplot(3,1,3)
histogram(Minv)
title('Invariant Mass (GeV/c^2)')
%% INVESTIGATE THE IMPACT OF IMPOSING MuQuality SELECTIONS CRITERIA (CUTS)
% ON THE DiMuon INVARIANT MASS DISTRIBUTION

%
% Histogram the DiMuon invariant mass for MuQuality >=0, MuQuality >=1,
% MuQuality >=2 and MuQuality >=3. Arrange the four histograms in
% one Figure using subplots.
%
IM0=find(MuQuality==0);
IM1=find(MuQuality==1);
IM2=find(MuQuality==2);
IM3=find(MuQuality==3);
figure('Name','Invariant mass for Dimuon events seperated by quality')
subplot(2,2,1)
histogram(Minv(IM0),0:0.33:7)
title('Invariant Mass by Quality: Quality 0 (GeV/c^2)')
subplot(2,2,2)
histogram(Minv(IM1),0:0.33:7)
title('Quality 1 (GeV/c^2)')
subplot(2,2,3)
histogram(Minv(IM2),0:0.33:7)
title('Quality 2 (GeV/c^2)')
subplot(2,2,4)
histogram(Minv(IM3),0:0.33:7)
title('Quality 3 (GeV/c^2)')

%% FOR EACH EVENT TRANSFORM TO THE REST FRAME OF THE DiMuon SYSTEM
% 1. TO CALCULATE THE DiMuon INVARIANT MASS IN THE REST FRAMEb
% 2. TO INVESTIGATE THE RELATIONSHIPS BETWEEN THE COMPONENTS OF THE
% 4-MOMENTUM OF THE FIRST MUON, [E1,px1,py1,pz1], WITH THAT OF THE
% SECOND MUON, [E2,px2,py2,pz2], IN THE REST FRAME
%
Bsys=[P3sys(:,1)./Esys,P3sys(:,2)./Esys,P3sys(:,3)./Esys];
%
B=sqrt(Bsys(:,1).^2+Bsys(:,2).^2+Bsys(:,3).^2);
Gam=sqrt(1./(1-B.^2));

% For the Lorentz Transformation ( http://en.wikipedia.org/wiki/Lorentz_transformation )
% we need to calculate the following quantities for the DiMuon Rest Frame:
% The DiMuon system's velocity: Bsys=(Bx,By,Bz)=P3sys/Esys;
% B=sqrt(Bx^2+By^2+Bz^2) and Gam=sqrt(1/(1-B^2))
%

%% TRANSFORM TO THE DiMuon REST FRAME - LOOP OVER THE EVENTS

    for i=1:Nevents
    LorTrans=[Gam(i) -Gam(i)*B(i) 0 0; -Gam(i)*B(i) Gam(i) 0 0; 0 0 1 0; 0 0 0 1];
    P4Mu1RF(i,:)=P4Mu1(i,:)*LorTrans;
    P4Mu2RF(i,:)=P4Mu2(i,:)*LorTrans;
    end
EsysRF=P4Mu1RF(:,1)+P4Mu2RF(:,1);
P3sysRF=[P4Mu1RF(:,2)+P4Mu2RF(:,2)  ,  P4Mu1RF(:,3)+P4Mu2RF(:,3) ,  P4Mu1RF(:,4)+P4Mu2RF(:,4)];
MinvRF=sqrt(EsysRF.^2-sum(P3sysRF.^2,2));
   

%% Histograms    
%
% Diagnostic Plots: Histogram all 10 input variables
% in the DiMuon Rest Frame in one Figure using subplots
%

figure('Name','Initial vars adjusted for Rest Frame')
subplot(3,4,1)
histogram(EvtNum)
title('Event Numbers')

% subplot(3,4,2)

subplot(3,4,3)
histogram(MuQuality,4)
title('Quality of Muon Data')


% subplot(3,4,4)

subplot(3,4,5)
histogram(P4Mu1RF(:,1))
title('E1 Energy Rest Frame (GeV)')


subplot(3,4,6)
histogram(P4Mu1RF(:,2))
title('X-component for P1 in RF (GeV/c)')

subplot(3,4,7)
histogram(P4Mu1RF(:,3))
title('Y-component for P1 in RF (GeV/c)')

subplot(3,4,8)
histogram(P4Mu1RF(:,4))
title('Z-component for P1 in RF (GeV/c)')

subplot(3,4,9)
histogram(P4Mu1RF(:,1))
title('E2 Energy in RF (GeV)')

subplot(3,4,10)
histogram(P4Mu1RF(:,2))
title('X-component for P2 in RF (GeV/c)')

subplot(3,4,11)
histogram(P4Mu1RF(:,3))
title('Y-component for P2 in RF (GeV/c)')

subplot(3,4,12)
histogram(P4Mu2RF(:,4))
title('Z-component for P2 in RF (GeV/c)')


%
% Histogram Esys, Psys (the magnitude of the 3-momentum)and the invariant
% mass of the DiMuon System in the DiMuon Rest Frame in one Figure using
% subplots
PsysRF=sqrt(P3sysRF(:,1).^2+P3sysRF(:,2).^2+P3sysRF(:,3).^2);
figure
subplot(3,1,1)
histogram(EsysRF)
title('Energy of the System Rest Frame (GeV)')
subplot(3,1,2)
histogram(PsysRF)
title('Magnitude of the 3-momentum in Rest Frame (GeV/c)')
subplot(3,1,3)
histogram(MinvRF)
title('Invariant Mass in Rest Frame (GeV/c^2)')

%
% Histogram the DiMuon Invariant mass in the DiMuon Rest Frame
% for MuQuality >=0, MuQuality >=1, MuQuality >=2 and MuQuality >=3.
% Arrange the four histograms in one Figure using subplots
%

figure
subplot(2,2,1)
histogram(MinvRF(IM0),0:0.33:7)
title('Invariant Mass by Quality in Rest Frame: Quality 0 (GeV/c^2)')
subplot(2,2,2)
histogram(MinvRF(IM1),0:0.33:7)
title('Quality 1 in RF (GeV/c^2)')
subplot(2,2,3)
histogram(MinvRF(IM2),0:0.33:7)
title('Quality 2 in RF (GeV/c^2)')
subplot(2,2,4)
histogram(MinvRF(IM3),0:0.33:7)
title('Quality 3 in RF (GeV/c^2)')
