function NStruct = NoiseTotalPoint(m1,h1,lat,lon,frq11,env,fpath)
%Function calculates total noise parameters such as median noise figure, 
%upper and lower deciles. Also it calculates the same parameters for the 
%noise components: atmospheric, man-made and galactic. 

%Inputs:
%m1 - month
%h1 - time in hours, LT
%lat - latitude, scalar, grad
%lon - longitude, scalar, grad
%frq11 - frequency, MHz
%env - man-made noise environment(1 - City, 2 - Residental, 3 - Rural, 4 - Quiet rural)
%fpath - path to the directory with Noise.mat file
%Outputs:
%NStruct - structure of noise parameters
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
        mustBeInRange(frq11,0.01,30)}
    env (1,1) {mustBeNumeric,mustBeMember(env,1:4)};
    fpath {mustBeText}
end

%Boundary points and centers of time blocks
h0 = [0,2:4:22,24];

%Structure initialization
NStruct.FmA = []; NStruct.DuA = []; NStruct.DlA = [];
NStruct.FmM = []; NStruct.DuM = []; NStruct.DlM = [];
NStruct.FmG = []; NStruct.DuG = []; NStruct.DlG = [];
NStruct.FmT = []; NStruct.DuT = []; NStruct.DlT = [];

%Median noise figure for the atmospheric component
FmA = NoiseFamPoint(m1,h0,lat,lon,frq11,fpath);
%Statistical parameters for the atmospheric component
[DuA,DlA,sigmaA_u,sigmaA_l] = DParam(m1,h0,lat,frq11,fpath);

%Parameters of man-made and galactic components
[FmM,DlM,DuM,FmG,DuG,DlG] = NoiseMG(env,frq11);

%Galactic noise standard deviation
sigmaG = 1.56;
%Man-made noise standard deviation
sigmaM_l = DlM/1.282;
sigmaM_u = DuM/1.282;

%c index
c = 10/log(10);

%alpha indexes
alphaA_u = exp(FmA/c+(sigmaA_u.^2)/2/(c^2));
alphaM_u = exp(FmM/c+(sigmaM_u^2)/2/(c^2));
alphaG = exp(FmG/c+(sigmaG^2)/2/(c^2)); %The same above and below median
alphaT_u = alphaA_u + alphaM_u + alphaG;

alphaA_l = exp(FmA/c+(sigmaA_l.^2)/2/(c^2));
alphaM_l = exp(FmM/c+(sigmaM_l^2)/2/(c^2));
alphaT_l = alphaA_l + alphaM_l + alphaG;

%beta indexes
betaA_u = (alphaA_u.^2).*(exp((sigmaA_u/c).^2)-1);
betaM_u = (alphaM_u^2)*(exp((sigmaM_u/c).^2)-1);
betaG = (alphaG^2)*(exp((sigmaG/c)^2)-1); %The same above and below median
betaT_u = betaA_u + betaM_u + betaG;

betaA_l = (alphaA_l.^2).*(exp((sigmaA_l/c).^2)-1);
betaM_l = (alphaM_l^2)*(exp((sigmaM_l/c).^2)-1);
betaT_l = betaA_l + betaM_l + betaG;

%gamma index
gammaT = exp(FmA/c)+exp(FmM/c)+exp(FmG/c);

%Total noise standard deviation, above median
if any(DuA>12,'all')
    sigmaT_u = zeros(size(DuA));
    sigmaT_u(DuA>12) = c*sqrt(2*log(alphaT_u(DuA>12)./gammaT(DuA>12)));
    sigmaT_u(DuA<=12) = c*sqrt(log(1 + betaT_u(DuA<=12)./(alphaT_u(DuA<=12).^2)));
else
    sigmaT_u = c*sqrt(log(1 + betaT_u./(alphaT_u.^2)));
end

%Total noise standard deviation, below median
if any(DlA>12,'all')
    sigmaT_l = zeros(size(DlA));
    sigmaT_l(DlA>12) = c*sqrt(2*log(alphaT_l(DlA>12)./gammaT(DlA>12)));
    sigmaT_l(DlA<=12) = c*sqrt(log(1 + betaT_l(DlA<=12)./(alphaT_l(DlA<=12).^2)));
else
    sigmaT_l = c*sqrt(log(1 + betaT_l./(alphaT_l.^2)));
end

%Total noise upper and lower deciles
DlT = sigmaT_l*1.282;
DlT1 = interp1(h0,DlT,h1,'spline');
DuT = sigmaT_u*1.282;
DuT1 = interp1(h0,DuT,h1,'spline');

FamT_l =  c*(log(alphaT_l) - (sigmaT_l.^2)/(2*c^2));
FamT_u =  c*(log(alphaT_u) - (sigmaT_u.^2)/(2*c^2));

FmT = max(FamT_l,FamT_u);
FmT1 = interp1(h0,FmT,h1,'spline');

%write outputs
NStruct.FmT = FmT1; %Worst-case noise
NStruct.DlT = DlT1;
NStruct.DuT = DuT1;
NStruct.FmA = FmA;
NStruct.DlA = DlA;
NStruct.DuA = DuA;
NStruct.FmM = FmM;
NStruct.DlM = DlM;
NStruct.DuM = DuM;
NStruct.FmG = FmG;
NStruct.DlG = DlG;
NStruct.DuG = DuG;