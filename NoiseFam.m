function Fam = NoiseFam(m1,h1,lat,lon,frq11,fpath,IntSz,options)
%Function calculates atmospheric noise figure approximated to the chosen 
%frequency 0.003-30. Threre is a linear extrapolation in 0.003-0.01 MHz
%frequency band. Also function performs equatorial interpolation to smooth
%transequatorial discontinuity

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

%Output:
%Fam - median atmospheric noise figure in dB above kT0b

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

%Проверка наличия файла
if ~(exist(fullfile(fpath,'Noise.mat'),'file') == 2)
    error("Can't find 'Noise.mat' in "+fpath)
end

%Загружаем интерполяторы
load(fullfile(fpath,'Noise.mat'),'F1','F2')
%Координатная сетка
[lat1,lon1] = ndgrid(lat,lon);
if options.time == "UT"
    %Пересчет времени UT в LT
    h1 = mod(h1 + lon1/15,24);
else
    %Оставляем время LT
    h1 = h1*ones(size(lon1));
end

%Изменяем размер m1
m1 = m1*ones(size(lat1));
%Используем первый интерполятор
%Получаем коэффициент шума на частоте 1 МГц
Fa3 = F1(m1,h1,lat1,lon1);

if options.interpolate == true %если интерполяция включена

    if max(lat) < 0 %Весь район в южном полушарии
        %Интерполяция не нужна
        %Месяц m1 со сдвигом на полгода
        Fam = F2(mod(m1+6,12),h1,Fa3,log10(frq11*ones(size(Fa3))));
    elseif min(lat) > 0 %Весь район в северном полушарии
        %Интерполяция не нужна
        %Месяц m1 без сдвига на полгода
        Fam = F2(m1,h1,Fa3,log10(frq11*ones(size(Fa3))));
    else %Район захватывает экватор
        
        %Приведение к заданной частоте (южное полушарие, сдвиг на 6 мес.)
        Fam1 = F2(mod(m1+6,12),h1,Fa3,log10(frq11*ones(size(Fa3))));
        %Приведение к заданной частоте (северное полушарие)
        Fam2 = F2(m1,h1,Fa3,log10(frq11*ones(size(Fa3))));
        
        %Индексы вектора широт
        latIdx = 1:numel(lat);
        %Индексы широт южниой части южного полушария (без интерполяции)
        latIdxS = latIdx(lat<0-IntSz/2);
        %Северная часть северного полушария (без интерполяции)
        latIdxN = latIdx(lat>0+IntSz/2);
        %Индексы полосы интерполяции
        latIdxI = latIdx(lat>=0-IntSz/2 & lat<=0+IntSz/2);
        
        %Весовые коэффициенты
        weights = (0+1/(numel(latIdxI)+1):1/(numel(latIdxI)+1):1-1/(numel(latIdxI)+1))';
        
        %Коэффициент шума с интерполяцией
        Fam = [Fam1(latIdxS,:);...
            Fam1(latIdxI,:).*flip(weights) + Fam2(latIdxI,:).*weights;...
            Fam2(latIdxN,:)];
    end
else %если интерполяция отключена
    if max(lat) < 0 %Весь район в южном полушарии
        %Интерполяция не нужна
        %Месяц m1 со сдвигом на полгода
        Fam = F2(mod(m1+6,12),h1,Fa3,log10(frq11*ones(size(Fa3))));
    elseif min(lat) > 0 %Весь район в северном полушарии
        %Интерполяция не нужна
        %Месяц m1 без сдвига на полгода
        Fam = F2(m1,h1,Fa3,log10(frq11*ones(size(Fa3))));
    else %Район захватывает экватор
        %Приведение к заданной частоте (южное полушарие, сдвиг на 6 мес.)
        Fam1 = F2(mod(m1+6,12),h1,Fa3,log10(frq11*ones(size(Fa3))));
        %Приведение к заданной частоте (северное полушарие)
        Fam2 = F2(m1,h1,Fa3,log10(frq11*ones(size(Fa3))));

        %Индексы вектора широт
        latIdx = 1:numel(lat);
        %Индексы широт южниой части южного полушария (без интерполяции)
        latIdxS = latIdx(lat<0);
        %Северная часть северного полушария (без интерполяции)
        latIdxN = latIdx(lat>=0);

        %Коэффициент шума без интерполяции
        Fam = [Fam1(latIdxS,:);Fam2(latIdxN,:)];
    end

end