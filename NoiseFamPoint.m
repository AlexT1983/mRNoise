function Fam = NoiseFamPoint(m,h,lat,lon,frq11,fpath)
%Function calculates median atmospheric noise figure for a given geographic
%point, month and time. Between middle timepoints of time blocks 2:4:22
%function performs spline interpolation. Calculations avaolable only fo LT
%time.

%Output:
%Fam - matrix where rows correspond to the elements of m argument,
% columns - to the elements of h argument

%m - month, vector
%h - local time, vector
%lat - latitude, scalar, grad
%lon - longitude, scalar, grad
%frq11 - frequency, MHz
%fpath - path to the directory with Noise.mat file

arguments
    m {mustBeNumeric,mustBeReal,mustBeVector,...
        mustBeInRange(m,0,12)}
    h {mustBeNumeric,mustBeReal,mustBeVector,...
        mustBeInRange(h,0,24)}
    lat (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(lat,-90,90)}
    lon (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(lon,-180,180)}
    frq11 (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(frq11,0.01,30)}
    fpath {mustBeText}
end

%Check the Noise.mat file existance
if ~(exist(fullfile(fpath,'Noise.mat'),'file') == 2)
    error("Can't find 'Noise.mat' in "+fpath)
end

%Calculations perfomed for the local time
%Middle points of time blocks
h0 = 2:4:22;
h1 = [0,h0,24];

%Load interpolators
load(fullfile(fpath,'Noise.mat'),'F1','F2')

%Fullgrid month/hour
[m2,h2] = ndgrid(m,h1);

%Change size of the lat and lon arguments
lat2 = lat*ones(size(m2));
lon2 = lon*ones(size(m2));

%Use first interpolator
Fa3 = F1(m2,h2,lat2,lon2);

if lat < 0 %Geographic point is in the southern hemisphere
    %Shift m2 for 6 months
    Fam1 = F2(mod(m2+6,12),h2,Fa3,log10(frq11*ones(size(Fa3))));
else %Geographic point is in the northern hemisphere
    %Don't shift m2 for 6 months
    Fam1 = F2(m2,h2,Fa3,log10(frq11*ones(size(Fa3))));
end

Fam = zeros(numel(m),numel(h));
%Spline interpolation of the timepoints
for i = 1:numel(m)
    Fam(i,:) = interp1(h1,Fam1(i,:),h,"spline");
end

end