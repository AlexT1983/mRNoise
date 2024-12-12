function [P,E] = NoisePwrPoint(m1,h1,lat,lon,frq11,env,...
    bandwidth,p,fpath)
%Function performs calculation of noise power and noise E-field

%Inputs:
%m1 - month
%h1 - time in hours, LT
%lat - latitude, grad
%lon - longitude, grad
%frq11 - frequency for Fam calculation, MHz
%env - man-made noise environment(1 - City, 2 - Residental, 3 - Rural, 4 - Quiet rural)
%bandwidth - frequency bandwidth, Hz
%p - probability that noise will not exceed given value
%Outputs:
%P - noise power in dB abovw 1 Watt
%E - noise E-field in dB above 1 uVolt/meter
%assumed short receiving monopole antenna on the perfect ground

arguments
    m1 {mustBeNumeric,mustBeReal,mustBeVector,...
        mustBeInRange(m1,0,12)}
    h1 {mustBeNumeric,mustBeReal,mustBeVector,...
        mustBeInRange(h1,0,24)}
    lat (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(lat,-90,90)}
    lon (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(lon,-180,180)}
    frq11 (1,1) {mustBeNumeric,mustBeReal,...
        mustBeInRange(frq11,0.003,30)}
    env (1,1) {mustBeNumeric,mustBeMember(env,1:4)};
    bandwidth (1,1) {mustBeNumeric,mustBeReal,...
        mustBePositive}
    p (1,1) {mustBeNumeric,mustBePositive,...
        mustBeInRange(p,0,1,'exclusive')}
    fpath {mustBeText}
end

NStruct = NoiseTotalPoint(m1,h1,lat,lon,frq11,env,fpath);

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