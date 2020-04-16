%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to test CSP formation for EOM
% Use simple straight ray homogenoues travel times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

no=1;
dt=4e-3;
dx=12.5;
Nr=640;
Ns=160;
vel=1480;
Lsec=6.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locS(:,1)=[0:4*dx:4*dx*Ns];
locS(:,2)=0;
locR(:,1)=[0:dx:dx*Nr];
locR(:,2)=0;
locP=[dx*Nr/2 2000];
wave=ricker(1/dt/2/10,dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(no);clf;
set(gcf, 'Color', [1 1 1]);
get(0,'ScreenSize');
set(0,'Units','centimeters');
set(gcf,'units','centimeters');
%[left bottom width height] A4 21 x 29.5 cm
pos=[72 30 24 18];
set(gcf,'position',pos);
set(gcf,'paperpositionmode','auto');


plot(locR(:,1),locR(:,2),'rv');
hold on;
plot(locS(:,1),locS(:,2),'k*');
plot(locP(:,1),locP(:,2),'bsq');
axis equal
hold off;
set(gca,'ydir','reverse');
box on;
%axis([0 (Nr+1)*12.5 0 max(tt(1:L))]);
%xlabel('Offset (m)','fontsize',14);
%ylabel('Time (s)','fontsize',14);
%title('Shot Gather','fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=Lsec/dt;
tt(1:L)=[1:L]*dt-dt;

for s=1:Ns
for r=1:Nr

locP(1,1)=(locS(s,1)+locR(r,1))/2; % making constant hoirzon reflector
srcd = sqrt( (locS(s,1)-locP(1,1) )^2 + (locS(s,2)-locP(1,2) )^2 );
recd = sqrt( (locR(r,1)-locP(1,1) )^2 + (locR(r,2)-locP(1,2) )^2 );
G(1:L)=0;
ti=round( (srcd+recd)/vel/dt );
G(ti)=1/(srcd+recd);
aa=conv(G,wave);
trace(s,r,1:L)=aa(1:L);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save test_layer_matlab trace
%load test_layer_matlab.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1;
il(1:Ns*Nr,1)=0;
xl(1:Ns*Nr,1)=0;
X(1:Ns*Nr,1)=0;
Y(1:Ns*Nr,1)=0;
cdp(1:Ns*Nr,1)=0;
seismic(1:Ns*Nr,1:L)=0;

for s=1:Ns
for r=1:Nr
%il(s,r)=locR(r,1);
%xl(s,r)=locR(r,2);
%X(s,r)=(locR(r,1)+locS(s,1))/2;
%Y(s,r)=(locR(r,2)+locS(s,2))/2;
sx(k,1)=locS(s,1);
sy(k,1)=locS(s,2);
gx(k,1)=locR(r,1);
gy(k,1)=locR(r,2);

il(k,1)=locR(r,1);
xl(k,1)=locR(r,2);
X(k,1)=(locR(r,1)+locS(s,1))/2;
Y(k,1)=(locR(r,2)+locS(s,2))/2;

%2d eqns
cdp(k,1)=(locR(r,1)+locS(s,1))/2/6.25;
off(k,1)=abs( locR(r,1)-locS(s,1) );

seismic(k,1:L)=trace(s,r,:);
k=k+1;
end
end

%WriteSegy('test_layer_matlab.segy',seismic','dt',0.004,'Inline3D',il','Crossline3D',xl','cdpX',X','cdpY',Y');

%write_outsegy('layer_ibm.segy',seismic','dt',0.004,'Inline3D',il','Crossline3D',xl', ...
%              'cdpX',X','cdpY',Y','SX',sx','SY',sy','GX',gx','GY',gy','CDP',cdp','offset',off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(no+1);clf;
set(gcf, 'Color', [1 1 1]);
get(0,'ScreenSize');
set(0,'Units','centimeters');
set(gcf,'units','centimeters');
%[left bottom width height] A4 21 x 29.5 cm
pos=[72 3 24 18];
set(gcf,'position',pos);
set(gcf,'paperpositionmode','auto');

s=10;

mm=max(max(max(trace)));
for r=1:10:Nr
tr(:,1)=trace(s,r,:);
plot(10*tr/mm+r,tt,'b');
hold on;
set(gca,'ydir','reverse');
box on;
end
%axis([0 (Nr+1)*12.5 0 max(tt(1:L))]);
%xlabel('Offset (m)','fontsize',14);
%ylabel('Time (s)','fontsize',14);
title('Shot Gather','fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

clear all
%model thickness(km) density(g/cm^3) vp(km/s)  Qp  vs(km/s)  Qs   
m1 = [0.25    1.8   2.0   1200    1.5   1000
      0.50    2.0   2.5   1200    2.0   1000];

vs(1)=1.5;
%%%%%%%%%%%%%%%%

geo=[0.0:0.25:2.0];
nstat=length(geo);
z_Source = .005;

%%%%%%%%%%%%%%%%

Tmax=6;                         % Desired record length in seconds
dt=0.004;                       % Time sample rate (s)
fmin=30;                        % Frequency band for ref modeling
fmax=40;                        % Frequency band for ref modeling
smin=0;                         % Minimum slowness
ds=0.01/fmax*max(geo);          % Slowness increment
smax=round(1.2/vs(1)/ds)*ds;    % Maximum slowness
tau=15;                         % for complex freq to supress wrap around
per_s=100*(1-1/1.2);            % percentual of the slowness integral that is windowed
per_f=40;                       % parameter for Hanning window  percentual of the FREQUECY that is windowed
direct =1;                      % to compute directwave = 1, else = 0
multiple=0  ;                   % ro compute multiples = 1, else = 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters=[Tmax           %T_max:   maximum modeling time (sec);
dt           %delta_T: time interval (msec);
smin         %u1:      initial slowness to compute (sec/km); 
smax         %u2:      last slowness to comute (sec/km);
fmin         %f1:      initial frequency to model (Hz);
fmax         %f2:      last frequency to model (Hz);
tau            %tau:     wrap around atenuation;

z_Source      % Depth source (km);
per_s           %perc_U   percentual of the slowness integral that is windowed
per_f           %perc_F   percentual of the FREQUECY that is windowed
0            % F1     '1'  the single force cartesian components (Equation (62) 
0            % F2      '1' from Muller, 1985)
1            % F3    '10'
direct            % to compute directwave = 1, else = 0
multiple            % ro compute multiples = 1, else = 0
0.0001];          % delta_U: slowness interval (0.0006 for offset 1.0 km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[uR,uZ,t] = reflectivity(m1,parameters,geo);

figure(3);clf
for i=1:nstat
plot(t',i+uZ(:,i)/max(uZ(:,i)));
hold on;
end
hold off

%figure;
%plotseis(1,uR2,t,geo);ylabel('Time (second)');xlabel('Offset (km)');title('radial component')
%figure;
%plotseis(1,uZ2,t,geo);ylabel('Time (second)');xlabel('Offset (km)');title('vertical component')

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all

    nlay = 2;  % Number of layers including halfspace:
    thick = [ 0.5   0.5 ]; %   0.5     0.5     0.0];
    alpha = [ 1.5   2.0 ]; %   2.5     3.0     3.6];
    beta  = [ .78   1.2 ]; %   1.3     1.7     1.9];
    den   = [ 1.9   2.0 ]; %  2.1     2.2     2.3];

    direct = 0;         % direct = 0 ==> include direct: 1 ==> delete direct:
    tmax = 2.5;         % Trace length in seconds:
    f0 = 30.0;          % Predominant frequency (Hz)of Gabor wavelet:
    gam = 4.0;          % Waist parameter of Gabor wavelet: Hardwired:

    rdist1 = 0.0;       % Initial offset in km:
    rstep  = 0.05;      % Distance between adjacent receivers in km:
    nrdist = 101;       % Number of receivers:

    nsdepth = 5;        % Source depth in grid points Receiver depth "nrdepth" hardwired at Surface Receivers:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_in = [ thick alpha beta den ]';
offset_in = [rdist1 rstep nrdist nsdepth]';
rdoff = rdist1 + rstep*(( 1 : nrdist ) - 1); % From offset_in

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate distance in m or km to reflecting pseudo-boundary:
%
    aar = 2.0*max(rdoff);        % aar = 3.0*max(rdoff)
                                 % THIS DOESN'T HAPPEN TOO OFTEN:
                                 % Only for long surface arrays
                                 % 2.0*max(rdoff) is usually enough
                                 % to prevent spurious reflections from
                                 % the reflecting pseudo-boundary at
                                 % r=a: It is instructive to run this
                                 % model with aar = 2.0*max(rdoff):
                                 % As alpha(1)=2.0, the total number of
                                 % roots required to be considered in the
                                 % Wavelength_Period Domain is thus
                                 % 2.0*(3.0*max(rdoff)):
    ppwl = 20;                   % Points Per Wavelet --> 20 or 40:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main loop

    para_in = [tmax, f0, gam, aar nlay ppwl direct]';

    [vv,hh,dt,jtime,ntimp ] = reflectfd_ioff(model_in,offset_in,para_in);

    [mtime,noff] = size(vv);                  % Just checking:
    T2s   = 1.0/f0;                           % Periods to seconds:
    atime = dt*jtime*( 1 : mtime )*T2s;       % Time axis definition.
    atime = atime';

nstat=101;
figure(8);clf
for i=1:nstat
plot(i+vv(:,i)/max(vv(:,i)));
hold on;
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
