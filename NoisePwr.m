function [P,E] = NoisePwr(m1,h1,lat,lon,frq11,env,...
    bandwidth,p,fpath,IntSz,options)
%Function performs calculation of noise power and noise E-field

%Inputs:
%m1 - month
%h1 - time in hours, UT or LT
%lat - vector of latitudes, grad
%lon - vector of longitudes, grad
%frq11 - frequency for Fam calculation, MHz
%env - man-made noise environment(1 - City, 2 - Residental, 3 - Rural, 4 - Quiet rural)
%bandwidth - frequency bandwidth, Hz
%p - probability that noise will not exceed given value
%fpath - path to the directory with Noise.mat file
%IntSz - width of the equatorial interpolation zone, grad (optional argument)
%Options:
%'time','LT/'UT' - local time or universal coordinated time
%'interpolation',true/false - turn on/off interpolation
%Outputs
%P - noise power in dB abovw 1 Watt
%E - noise E-field in dB above 1 uVolt/meter
%assumed short receiving monopole antenna on the perfect ground

arguments
    m1 (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(m1,0,12)}
    h1 (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(h1,0,24)}
    lat {mustBeNumeric,mustBeReal,...
        mustBeVector,...
        mustBeInRange(lat,-90,90)}
    lon {mustBeNumeric,mustBeReal,...
        mustBeVector,...
        mustBeInRange(lon,-180,180)}
    frq11 (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(frq11,0.003,30)}
    env (1,1) {mustBeNumeric,mustBeMember(env,1:4)};
    bandwidth (1,1) {mustBeNumeric,mustBeReal,...
        mustBePositive}
    p (1,1) {mustBeNumeric,mustBePositive,...
        mustBeInRange(p,0,1,'exclusive')}
    fpath {mustBeText}
    IntSz (1,1) {mustBeNumeric,mustBePositive,...
        mustBeInRange(IntSz,10,50)} = 20;
    options.time {mustBeMember(options.time,{'UT','LT'})} = 'LT'
    options.interpolate (1,1) logical = false
end

NStruct = NoiseTotal(m1,h1,lat,lon,frq11,env,fpath,IntSz,...
    'time',options.time,'interpolate',options.interpolate);

if p >= 0.5
    E = NStruct.FmT + NStruct.DuT + 20*log10(frq11)...
        + 10*log10(bandwidth) - 95.5; %dB(uV/m)
    P = NStruct.FmT + NStruct.DuT + 10*log10(bandwidth) - 204; %dBW
elseif p == 0.5 %Use median noise figure
    E = NStruct.FmT + 20*log10(frq11)...
        + 10*log10(bandwidth) - 95.5; %dB(uV/m)
    P = NStruct.FmT + 10*log10(bandwidth) - 204; %dBW
else
    E = NStruct.FmT - NStruct.DlT + 20*log10(frq11)...
        + 10*log10(bandwidth) - 95.5; %dB(uV/m)
    P = NStruct.FmT - NStruct.DlT + 10*log10(bandwidth) - 204; %dBW
end