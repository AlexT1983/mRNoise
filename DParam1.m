function [Dua,Dla,sigma_u,sigma_l] = DParam1(m1,h1,lat,lon,frq11,...
    fpath,IntSz,options)
%Function calculates lower and upper deciles of atmospheric noise figure
%distribution and corresponding standard deviations

%Inputs:
%m1 - month
%h1 - time in hours, UT or LT
%lat - vector of latitudes, grad
%lon - vector of longitudes, grad
%frq11 - frequency for Fam calculation, MHz
%fpath - path to the directory with Noise.mat file
%IntSz - width of the equatorial interpolation zone, grad (optional argument)
%Options:
%'time','LT/'UT' - local time or universal coordinated time
%'interpolation',true/false - turn on/off interpolation
%Outputs:
%Dua - atmospheric noise upper decile, dB
%Dla - atmospheric noise lower decile, dB
%sigma_u - atmospheric noise standard deviation of the upper distribution, dB
%sigma_l - atmospheric noise standard deviation of the lower distribution, dB

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
    fpath {mustBeText}
    IntSz (1,1) {mustBeNumeric,mustBePositive,...
        mustBeInRange(IntSz,10,50)} = 20;
    options.time {mustBeMember(options.time,{'UT','LT'})} = 'LT'
    options.interpolate (1,1) logical = false
end

%Coordinates grid
[lat1,lon1] = ndgrid(lat,lon);
if options.time == "UT"
    %conver UT to LT
    h1 = mod(h1 + lon1/15,24);
else
    %Leave LT
    h1 = h1*ones(size(lon1));
end
%Change m1 size
m1 = m1*ones(size(lat1));
%find the logarithm of the frequency and change size
logfreq1 = log10(frq11)*ones(size(lon1));

%Make standard normal distribution
pd = makedist("Normal");

%Load interpolators
load(fullfile(fpath,'Noise.mat'),'FDuA','FDlA')

if options.interpolate == true %if interpolation enabled

    if max(lat) < 0 %Whole area in southern hemisphere
        %Interpolation isn't needed
        %Shift m1 for 6 month
        %Use interpolator for Du (p=0.9)
        Dua = abs(FDuA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_u = Dua/icdf(pd,0.9);

        %Use interpolator for Dl (p=0.1)
        Dla = abs(FDlA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_l = -Dla/icdf(pd,0.1);
    elseif min(lat) > 0 %Whole area in northern hemisphere
        %Interpolation isn't needed
        %Don't shift m1 for 6 month
        %Use interpolator for Du (p=0.9)
        Dua = abs(FDuA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_u = Dua/icdf(pd,0.9);

        %Use interpolator for Dl (p=0.1)
        Dla = abs(FDlA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_l = -Dla/icdf(pd,0.1);
    else %Equator is in the area
        %Southern hemisphere
        %Use interpolator for Du (p=0.9)
        Dua1 = abs(FDuA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_u1 = Dua1/icdf(pd,0.9);
        %Northern hemisphere
        %Use interpolator for Du (p=0.9)
        Dua2 = abs(FDuA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_u2 = Dua2/icdf(pd,0.9);

        %Southern hemisphere
        %Use interpolator for Dl (p=0.1)
        Dla1 = abs(FDlA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_l1 = Dla1/icdf(pd,0.1);
        %Northern hemisphere
        %Use interpolator for Dl (p=0.1)
        Dla2 = abs(FDlA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_l2 = Dla2/icdf(pd,0.1);


        %Latitude vector indexes
        latIdx = 1:numel(lat);
        %Latitude indexes for southern part of southern hemisphere (without interpolation)
        latIdxS = latIdx(lat<0-IntSz/2);
        %Nothern part of northern hemisphere (without interpolation)
        latIdxN = latIdx(lat>0+IntSz/2);
        %Interpolation area indexes
        latIdxI = latIdx(lat>=0-IntSz/2 & lat<=0+IntSz/2);
        
        %Weight indexes
        weights = (0+1/(numel(latIdxI)+1):1/(numel(latIdxI)+1):1-1/(numel(latIdxI)+1))';
        
        %Parameters Dua and Dla with interpolation
        Dua = [Dua1(latIdxS,:);...
            Dua1(latIdxI,:).*flip(weights) + Dua2(latIdxI,:).*weights;...
            Dua2(latIdxN,:)];
        Dla = [Dla1(latIdxS,:);...
            Dla1(latIdxI,:).*flip(weights) + Dla2(latIdxI,:).*weights;...
            Dla2(latIdxN,:)];
        %Parameters sigma_u and sigma_l with interpolation
        sigma_u = [sigma_u1(latIdxS,:);...
            sigma_u1(latIdxI,:).*flip(weights) + sigma_u2(latIdxI,:).*weights;...
            sigma_u2(latIdxN,:)];
        sigma_l = [sigma_l1(latIdxS,:);...
            sigma_l1(latIdxI,:).*flip(weights) + sigma_l2(latIdxI,:).*weights;...
            sigma_l2(latIdxN,:)];
    end
else %if interpolation is disabled
    if max(lat) < 0 %Whole area in southern hemisphere
        %Interpolation isn't needed
        %Shift m1 for 6 month
        %Use interpolator for Du (p=0.9)
        Dua = abs(FDuA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_u = Dua/icdf(pd,0.9);

        %Use interpolator for Dl (p=0.1)
        Dla = abs(FDlA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_l = -Dla/icdf(pd,0.1);
    elseif min(lat) > 0 %Whole area in southern hemisphere
        %Interpolation isn't needed
        %Don't shift m1 for 6 month
        %Use interpolator for Du (p=0.9)
        Dua = abs(FDuA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_u = Dua/icdf(pd,0.9);

        %Use interpolator for Dl (p=0.1)
        Dla = abs(FDlA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_l = -Dla/icdf(pd,0.1);
    else %Equator is in the area
        %Southern hemisphere
        %Use interpolator for Du (p=0.9)
        Dua1 = abs(FDuA(mod(m1+6,12),h1,logfreq1));
        %Find std for given probability
        sigma_u1 = Dua1/icdf(pd,0.9);
        %Northern hemisphere
        %Use interpolator for Du (p=0.9)
        Dua2 = abs(FDuA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_u2 = Dua2/icdf(pd,0.9);

        %Southern hemisphere
        %Use interpolator for Dl (p=0.1)
        Dla1 = abs(FDlA(mod(m1+6,12),h1,logfreq1));
        %НFind std for given probability
        sigma_l1 = Dla1/icdf(pd,0.1);
        %Northern hemisphere
        %Use interpolator for Dl (p=0.1)
        Dla2 = abs(FDlA(m1,h1,logfreq1));
        %Find std for given probability
        sigma_l2 = Dla2/icdf(pd,0.1);

        %Latitude vector indexes
        latIdx = 1:numel(lat);
        %Latitude indexes for southern part of southern hemisphere (without interpolation)
        latIdxS = latIdx(lat<0);
        %Nothern part of northern hemisphere (without interpolation)
        latIdxN = latIdx(lat>=0);

        %Parameters Dua и Dla without interpolation
        Dua = [Dua1(latIdxS,:); Dua2(latIdxN,:)];
        Dla = [Dla1(latIdxS,:); Dla2(latIdxN,:)];
        %Parameter sigma without interpolation
        sigma_u = [sigma_u1(latIdxS,:); sigma_u2(latIdxN,:)];
        sigma_l = [sigma_l1(latIdxS,:); sigma_l2(latIdxN,:)];
    end
end

end
