function [FmM,DlM,DuM,FmG,DuG,DlG] = NoiseMG(env,freq)

%Function calculates man-made and galactic noise parameters due to
%ITU-R P.372-17 Recommendation

%Outputs:
%FmM - median noise figure of the man-made noise
%DlM - lower decile of the man-made noise
%DuM - upper decile of the man-made noise
%FmG - median noise figure of the galactic noise
%DlG - lower decile of the galactic noise
%DuG - upper decile of the galactic noise

%Inputs:
%env - man-made noise environment (1 - city, 2 - residental,
%       3 - rural, 4 - quiet rural)
%freq - frequency, MHz

arguments
    env (1,1) {mustBeNumeric,mustBeMember(env,1:5)};
    freq {mustBeNumeric,mustBeVector,...
        mustBeInRange(freq,0.01,30)}
end

c = [76.8; 72.5; 67.2; 53.6; 52.0]; 
d = [27.7; 27.7; 27.7; 28.6; 23.0];

switch env
    case 1 %City
        FmM = c(1) - d(1)*log10(freq);
        DuM = 11.0;
        DlM = 6.7;
    case 2 %Residential
        FmM = c(2) - d(2)*log10(freq);
        DuM = 10.6;
        DlM = 5.3;
    case 3 %Rural
        FmM = c(3) - d(3)*log10(freq);
        DuM = 9.2;
        DlM = 4.6;
    case 4 %Quiet rural
        FmM = c(4) - d(4)*log10(freq);
        DuM = 9.2;
        DlM = 4.6;
end

FmG = c(5) - d(5)*log10(freq);
DlG = 2.0;
DuG = 2.0;