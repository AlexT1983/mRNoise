function [DuA,DlA,sigma_u,sigma_l] = DParam(m,h,lat,frq11,fpath)
%Function calculates lower and upper deciles of atmospheric noise figure
%distribution and corresponding standard deviations. Also it performs
%spline interpolation for the time points

%m - month, vector
%h - local time, vector
%lat - latitude, scalar, grad
%frq11 - frequency, MHz
%fpath - path to the directory with Noise.mat file

arguments
    m {mustBeNumeric,mustBeReal,mustBeVector,...
        mustBeInRange(m,0,12)}
    h {mustBeNumeric,mustBeReal,mustBeVector,...
        mustBeInRange(h,0,24)}
    lat (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(lat,-90,90)}
    frq11 (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(frq11,0.01,30)}
    fpath {mustBeText}
end

%Load interpolators
load(fullfile(fpath,'Noise.mat'),'FDuA','FDlA')
%Стандартное нормальное распределение
pd = makedist("Normal");

%Centers of time blocks and boundaries
h0 = [0,2:4:22,24];
%Fullgrid
[m1,h2] = ndgrid(m,h0);

%find the logarithm of the frequency and change size
logfreq1 = log10(frq11)*ones(size(m1));

if lat<0 %southern hemisphere
    %shift m1 for 6 months
    m1 = mod(m1+6,12);
    %Use interpolators for upper and lower parts of noise distribution
    DuA1 = abs(FDuA(m1,h2,logfreq1));
    DlA1 = abs(FDlA(m1,h2,logfreq1));
    %Find std via deciles
    sigma_u1 = DuA1/icdf(pd,0.9);
    sigma_l1 = -DlA1/icdf(pd,0.1);
    
else %northern hemisphere
    %don't shift m1
    %Use interpolators for upper and lower parts of noise distribution
    DuA1 = abs(FDuA(m1,h2,logfreq1));
    DlA1 = abs(FDlA(m1,h2,logfreq1));
    %Find std via deciles
    sigma_u1 = DuA1/icdf(pd,0.9);
    sigma_l1 = -DlA1/icdf(pd,0.1);
    
end

%Preallocate variables
DuA = zeros(numel(m),numel(h));
DlA = zeros(numel(m),numel(h));
sigma_u = zeros(numel(m),numel(h));
sigma_l = zeros(numel(m),numel(h));

%Spline interpolation of the timepoints
for i = 1:numel(m)
    DuA(i,:) = interp1(h0,DuA1(i,:),h,"spline");
    DlA(i,:) = interp1(h0,DlA1(i,:),h,"spline");
    sigma_u(i,:) = interp1(h0,sigma_u1(i,:),h,"spline");
    sigma_l(i,:) = interp1(h0,sigma_l1(i,:),h,"spline");
end