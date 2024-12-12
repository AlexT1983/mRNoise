%% F1 gridded interpolant usage
%Use F1 gridded interpolant to calculate median atmospheric noise figure at
%1 MHz frequency in decibels above kT0b

%Input data
%path to Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise\Noise.mat';
%month
m = 1; %january
%hour
h = 2; %local time
%latitude
lat = -90:90;
%longitude
lon = -180:180;
%full grid format
[lat1,lon1] = ndgrid(lat,lon);
m1 = m*ones(size(lat1));
h1 = h*ones(size(lat1));

%load interpolator object
load(fpath,'F1')

%Calculate median noise figure at 1 MHz
Fam = F1(m1,h1,lat1,lon1);

%Show noise figure on the worldmap (Mapping toolbox is needed)
figure
worldmap('World');
tightmap;

%Show Fam
contourfm(lat1,lon1,Fam,'LineStyle','none');
%Show colorbar and label it
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = '\itF_a_m\rm (1 MHz), dB/(\itkT\rm_0\itb\rm)';
%Show coastline
geoshow('landareas.shp','FaceColor','none');
%% F2 gridded interpolant usage
%Use F2 gridded interpolant to calculate median atmospheric noise figure at
%desirable frequency in decibels above kT0b
%Frequency range is 0.003-30 MHz

%Input data
%path to Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise\Noise.mat';
%month
m = 1; %january
%hour
h = 2; %local time
%latitude
lat = -90:90;
%longitude
lon = -180:180;
%Frequency, MHz
frq = 5;

%full grid format
[lat1,lon1] = ndgrid(lat,lon);
m1 = m*ones(size(lat1));
h1 = h*ones(size(lat1));
frq1 = frq*ones(size(lat1));

%load interpolator objects
load(fpath,'F1','F2')

%First calculate median noise figure at 1 MHz
Fam1 = F1(m1,h1,lat1,lon1);

%When we approximate frequency in southern hemisphere we must shift 'm'
% by 6 month. It caused by seasonal differences in
%approximation functions
m2 = m1;
m2(lat1<0) = mod(m1(lat1<0)+6,12);
%In F2 input frequency is logarithmic
Fam2 = F2(m2,h1,Fam1,log10(frq1));

%Show frequency approximated noise figure on the worldmap 
%(Mapping toolbox is needed)
figure
worldmap('World');
tightmap;

%Show Fam2
contourfm(lat1,lon1,Fam2,'LineStyle','none');
%You can see equatorial discontinuity due to seasonal differences in
%approximation functions

%Show colorbar and label it
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = '\itF_a_m\rm (5 MHz), dB/(\itkT\rm_0\itb\rm)';
%Show coastline
geoshow('landareas.shp','FaceColor','none');
%% Other interpolators usage
%Use FDuA, FDlA, FsigmaFaA, FsigmaDuA, FV_d, Fsigma_V_d
%gridded interpolants to calculate statistical parameters of the
%atmospheric noise at the desirable frequency

%Input data
%path to Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise\Noise.mat';
%month
m = 7; %July
%hour
h = 8; %local time
%Frequency, MHz
frq = logspace(log10(0.01),log10(30));

%load interpolators
load(fpath,'FDuA','FDlA','FsigmaFaA','FsigmaDuA',...
    'FsigmaDlA','FV_d','Fsigma_V_d');

%Transform input data
frqLog = log10(frq);
m1 = m*ones(size(frqLog));
h1 = h*ones(size(frqLog));

%Show atmospheric noise statistics for 8:00 LT July
figure
semilogx(frq,FDuA(m1,h1,frqLog),...
    'Marker','.','Color','b')
hold on
semilogx(frq,FDlA(m1,h1,frqLog),...
    'Marker','.','Color','r')
semilogx(frq,FsigmaFaA(m1,h1,frqLog),...
    'Marker','.','Color','g')
semilogx(frq,FsigmaDuA(m1,h1,frqLog),...
    'Marker','.','Color','k')
semilogx(frq,FsigmaDlA(m1,h1,frqLog),...
    'Marker','.','Color','m')
semilogx(frq,FV_d(m1,h1,frqLog),...
    'Marker','.','Color','c')
semilogx(frq,Fsigma_V_d(m1,h1,frqLog),...
    'Marker','.','Color','y')
hold off

%Configure axes
xlabel('Frequency, MHz')
ylabel('dB')
xlim([0.01,30])
grid on

%Show legend
legend('\itD_u_A', '\itD_l_A', '\it\sigma_F_A', '\it\sigma_D_u_A', ...
    '\it\sigma_D_l_A', '\itV_d', '\it\sigma_V_d')
%% Functions for atmospheric noise
%month
m1 = 7; %July
%time, hour
h1 = 6;
%frequency, MHz
frq = 15;
%path to the directory with the Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise';
%Width of equatorial interpolation zone
IntSz = 30;
%latitude and longitude vectors
lat = -89:89;
lon = -180:180;
%1-p - probability of noise level exceedance
p = 0.2;

%Median noise figure at 15 MHz. Interpolation turned on, time - UTC
Fam = NoiseFam(m1,h1,lat,lon,frq,fpath,IntSz,'time','UT','interpolate',1);
%D-parameter. Fam+D will not be exceeded with the 'p' probability
D = NoiseD(m1,h1,lat,lon,frq,fpath,p,IntSz,'time','UT','interpolate',1);

%Visualize maps
figure;
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
%full grid format for latitude and longitude
[lat1,lon1] = ndgrid(lat,lon);

nexttile;%----------------------------first tile---------------------------
Map1 = worldmap('World');
setm(Map1,'FontSize',8)
tightmap;

%Show Fam
contourfm(lat1,lon1,Fam,'LineStyle','none');
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = '\itF_a_m\rm (15 MHz), dB/(\itkT\rm_0\itb\rm)';
%Show coastline
geoshow('landareas.shp','FaceColor','none');

nexttile;%----------------------------second tile--------------------------
Map1 = worldmap('World');
setm(Map1,'FontSize',8)
tightmap;

%Show D
contourfm(lat1,lon1,D,'LineStyle','none');
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = ['\itD\rm (p = ',num2str(p),'), dB'];
%Show coastline
geoshow('landareas.shp','FaceColor','none');

%% Diurnal variability of atmospheric noise parameters
%These fuctions are useful for diurnal variations, because they perform
%spline interpolation for time points. But in can only be performed for
%single geographic point
m = 2:3:12; %month of year
h = 0:0.5:24; %time vector
%frequency, MHz
frq = 15;
%path to the directory with the Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise';
%latitude and longitude of geografic point
lat = 40;
lon = -50;
%1-p - probability of noise level exceedance
p = 0.2;

%'Point' functions
%Median noise figure at 15 MHz
Fam = NoiseFamPoint(m,h,lat,lon,frq,fpath);
%Noise parameters. Fam+D will not be exceeded with the 'p' probability
[~,~,~,sigma_l] = DParam(m,h,lat,frq,fpath);

%Find quantile (D) for the p-probability of non-exceedance
%make standard normal distribution
pd = makedist("Normal");
%p < 0, then use sigma_l from the lower distribution part
D = sigma_l*icdf(pd,p);

plot(h,Fam,'-',h,Fam+D,'--')
grid on
xlabel('Local time, hour')
xlim([0 24])
xticks(0:2:24)
ylabel('\itF_a_m\rm (15 MHz), dB/(\itkT\rm_0\itb\rm)')
legend('\itF_a_m\rm, February', '\itF_a_m\rm, May', ...
    '\itF_a_m\rm, August', '\itF_a_m\rm, November', ...
    '\itF_a_m+D\rm, February', '\itF_a_m+D\rm, May', ...
    '\itF_a_m+D\rm, August', '\itF_a_m+D\rm, November',...
    'Location','best','NumColumns',2)
%% Atmospheric noise distribution at the geographic point
%Rio Grande, Argentina
%Latitude and longitude 
lat = -53.78;
lon = -67.70;
%Time, hour
h = 18;
%Month - November
%Seasonal shift +6 month for southern hemisphere
%(we shift month only for interpolators arguments,
%when we use functions, such as NoiseFam(...), month shift isn't necessary
m = mod(10+6,12); 
%Frequency
frq = 15; %15 МГц
%full path (with name) to the Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise\Noise.mat';

%Load interpolators
load(fpath,'F1','F2','FDuA','FDlA','FsigmaFaA','FsigmaDuA','FsigmaDlA')

Fa1MHz = F1(m,h,lat,lon);
FaOnFrq = F2(m,h,Fa1MHz,log10(frq));
Du = FDuA(m,h,log10(frq));
Dl = FDlA(m,h,log10(frq));
sigmaFaA = FsigmaFaA(m,h,log10(frq));
sigmaDuA = FsigmaDuA(m,h,log10(frq));
sigmaDlA = FsigmaDlA(m,h,log10(frq));

%Find PDF values for 0.1, 0.5 and 0.9 quantiles
%Standard normal distribution
pd = makedist('Normal');
%Find quantiles with iCDF
psi1 = icdf(pd,0.1);
psi2 = icdf(pd,0.5);
psi3 = icdf(pd,0.9);
%Upper distribution part
sigma_u = Du/icdf(pd,0.9);
pd1 = makedist('Normal',FaOnFrq,sigma_u);
%Lower distribution part
sigma_l = Dl/icdf(pd,0.9);
pd2 = makedist('Normal',FaOnFrq,sigma_l);

%Plot estimate of Fa at given frequency
x1 = FaOnFrq-2*Dl:(2*Dl)/50:FaOnFrq;
x2 = FaOnFrq:(2*Du)/50:FaOnFrq+2*Du;

figure
yyaxis right
hold on
%Plot PDF
area(x1,pdf(pd2,x1),'FaceAlpha',0.2,'FaceColor',"#D95319")
area(x2,pdf(pd1,x2),'FaceAlpha',0.2,'FaceColor',"#0072BD")
hold off
ylabel('Probability density')
yyaxis left
%Plot CDF
hold on
plot(x1,cdf(pd2,x1),'LineStyle','-','LineWidth',1.5,'Color',"#D95319")
plot(x2,cdf(pd1,x2),'LineStyle','-','LineWidth',1.5,'Color',"#0072BD")
%Lower and upper deciles
plot(FaOnFrq-Dl,cdf(pd2,FaOnFrq-Dl),...
    'LineStyle','none','Marker','o','LineWidth',1,'Color',"#D95319")
plot(FaOnFrq+Du,cdf(pd1,FaOnFrq+Du),...
    'LineStyle','none','Marker','o','LineWidth',1,'Color',"#0072BD")
hold off
ylabel('Cumulative probability')
xline(FaOnFrq,'-','\itF_a_m\rm','LabelVerticalAlignment','bottom',...
    'LineWidth',1)
xline(FaOnFrq+Du,'--','\itF_a_m + D_u\rm','LabelVerticalAlignment','middle')
xline(FaOnFrq-Dl,'--','\itF_a_m - D_l\rm',...
    'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
yline(0.9,'--','\itp\rm = 0.9')
yline(0.1,'--','\itp\rm = 0.1')
xlabel(['\itF_a\rm (',num2str(frq),' MHz), dB/(\itkT\rm_0\itb\rm)'])
grid off
title({'Atmospheric noise figure distribution';...
    sprintf('latitude: %5.2f, longitude: %5.2f',lat,lon);...
    sprintf('local time: %i (hour), month: %i',h,mod(m+6,12))})
%% Total noise example
m = 2; %February
h = 9; %9:00
frq = 20; %Frequency, 20 MHz
%Path to directory with the Noise.mat file 
fpath = 'F:\Computations\Matlab\VLF\RF noise';
%Equatorial interpolation zone width, grad
IntSz = 20;
%Latitude and longitude vectors
lat = -90:90;
lon = -180:180;
%full grid coordinates
[lat1,lon1] = ndgrid(lat,lon);
%Man-made noise category
%1 - City, 2 - Residential,
%3 - Rural, 4 - Quiet rural
env = 3;
%Options:
%'interpolate',true/false - turn on/off interpolation
%'time','UT'/'LT' - local or universal coordinated time

%with interpolation, universal time
NStruct1 = NoiseTotal(m,h,lat,lon,frq,env,fpath,...
    IntSz,'interpolate',true,'time','UT');
%without interpolation, local time
NStruct2 = NoiseTotal(m,h,lat,lon,frq,env,fpath,...
    IntSz,'interpolate',false,'time','LT');

% NStruct - noise parameters structure:
% NStruct.FmT - total median noise figure;
% NStruct.DlT - total noise figure lower decile;
% NStruct.DuT - total noise figure upper decile;
% NStruct.FmA - atmospheric median noise figure;
% NStruct.DlA - atmospheric noise figure lower decile;
% NStruct.DuA - atmospheric noise figure upper decile;
% NStruct.FmM - man-made median noise figure;
% NStruct.DlM - man-made noise figure lower decile;
% NStruct.DuM - man-made noise figure upper decile;
% NStruct.FmG - galactic median noise figure;
% NStruct.DlG - galactic noise figure lower decile;
% NStruct.DuG - galactic noise figure upper decile;

%Visualize maps
figure;
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile;%----------------------------first tile---------------------------
Map1 = worldmap('World');
setm(Map1,'FontSize',8)
tightmap;

%Show total noise from NStruct1
contourfm(lat1,lon1,NStruct1.FmT,'LineStyle','none');
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = ['\itF_m_T\rm (',num2str(frq),...
    ' MHz), dB/(\itkT\rm_0\itb\rm)'];
%Show coastline
geoshow('landareas.shp','FaceColor','none');
title({'Total noise figure';'interpolation enabled, time - UT'})


nexttile;%----------------------------second tile--------------------------
Map1 = worldmap('World');
setm(Map1,'FontSize',8)
tightmap;

%Show total noise from NStruct2
contourfm(lat1,lon1,NStruct2.FmT,'LineStyle','none');
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = ['\itF_m_T\rm (',num2str(frq),...
    ' MHz), dB/(\itkT\rm_0\itb\rm)'];
%Show coastline
geoshow('landareas.shp','FaceColor','none');
title({'Total noise figure';'interpolation disabled, time - LT'})
%% Noise power and E-field
m = 5; %May
h = 11.5; %11:30
frq = 0.5; %Frequency, 500 kHz
bandwidth = 1000; %Bandwidth, Hz;
%Path to directory with the Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise';
%Equatorial interpolation zone width, grad
IntSz = 25;
%Latitude and longitude vectors
lat = -90:90;
lon = -180:180;
%full grid coordinates
[lat1,lon1] = ndgrid(lat,lon);
%Probability, that noise will not exceed given value
p = 0.8;
%Man-made noise category
%1 - City, 2 - Residential,
%3 - Rural, 4 - Quiet rural
env = 2;
%Options:
%'interpolate',true/false - turn on/off interpolation
%'time','UT'/'LT' - local or universal coordinated time

%Use noisePwr function (assumed short monopole on the perfect ground as
%receiving antenna)
%P - noise power, dBW
%E - noise E-field, dB(uV/m)
[P,E] = NoisePwr(m,h,lat,lon,frq,env,bandwidth,p,fpath,...
    IntSz,'interpolate',true,'time','UT');

%Visualize maps
figure;
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile;%----------------------------first tile---------------------------
Map1 = worldmap('World');
setm(Map1,'FontSize',8)
tightmap;

%Show power
contourfm(lat1,lon1,P,'LineStyle','none');
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = ['\itP\rm (',num2str(frq),' MHz), dBW'];
%Show coastline
geoshow('landareas.shp','FaceColor','none');
title('Noise power')


nexttile;%----------------------------second tile--------------------------
Map1 = worldmap('World');
setm(Map1,'FontSize',8)
tightmap;

%Show E-field
contourfm(lat1,lon1,E,'LineStyle','none');
C = contourcbar;
C.Location = 'southoutside';
C.Label.String = ['\itE\rm (',num2str(frq),' MHz), dB(\muV/m)'];
%Show coastline
geoshow('landareas.shp','FaceColor','none');
title('Noise E-field')
%% Total noise diurnal variation. Power and E-field
%centers of seasons
m1 = 6; %summer
%time vector
h1 = 0:0.5:24;
%frequency, MHz
frq = 3.7;
%path to the directory with the Noise.mat file
fpath = 'F:\Computations\Matlab\VLF\RF noise';
%latitude and longitude of geografic point
lat = 40;
lon = -50;
%Man-made noise category
%1 - City, 2 - Residential,
%3 - Rural, 4 - Quiet rural
env = 3;
bandwidth = 500; %Bandwidth, Hz;
%Probability, that noise will not exceed given value
p = 0.8;

%Structure of noise parameters
NStruct = NoiseTotalPoint(m1,h1,lat,lon,frq,env,fpath);
%Find std from deciles
pd = makedist('Normal');
%Std of lower and upper parts of total noise figure distribution
sigmaT_l = -NStruct.DlT/icdf(pd,0.1);
sigmaT_u = NStruct.DuT/icdf(pd,0.9);

%Noise power and E-field
[P,E] = NoisePwrPoint(m1,h1,lat,lon,frq,env,...
    bandwidth,p,fpath);

%Plot variables
%Noise figure and std
figure
yyaxis left
plot(h1,NStruct.FmT,'-')
yyaxis right
plot(h1,sigmaT_l,'--')
hold on
plot(h1,sigmaT_u,'-.')
xlim([0 24])
legend(['\itF_m_T\rm (',num2str(frq),' MHz), dB/(\itkT\rm_0\itb\rm)'],...
    '\it\sigma_u_T\rm, dB','\it\sigma_l_T\rm, dB', 'Location','best')
grid on

%power and E-field
figure
yyaxis left
plot(h1,P)
yyaxis right
plot(h1,E)
xlim([0 24])
legend(['\itP\rm (',num2str(frq),' MHz), dBW'],...
    ['\itE\rm (',num2str(frq),' MHz), dB(\muV/m)'], 'Location','best')
grid on