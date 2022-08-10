% Warmest pixel composite
%
% OPeNDAP Science Team
% Copyright 2008
% $Version 2.0.0$

IQual = input('Input the quality level (>=) to use to mask the data: 0-7 (0 means do not apply quality mask) (0) ');
if isempty(IQual)
    IQual = 0;
end
IYear = input('Input year of sequence to composite: ');
IStart = input('Input first day of sequence to composite: ');
IEnd = input('Input last day of sequence to composite: ');
ImagePlane = input('Input the image plane in which to plot night composite - daytime composite and day-night composites will be plotted in the subsequent image planes (1): ');
if isempty(ImagePlane)
    ImagePlane = 1;
end
MaxTemperature = input('Input the maximum temperature for plots (30): ');
if isempty(MaxTemperature)
    MaxTemperature = 30;
end

% Get number of zeros before day number in filename.
if IStart >= 100
    Filler = '';
else if IStart >= 10
        Filler = '0';
    else
        Filler = '00';
    end
end

% Mask data
if IQual > 0
    eval(['nnd1 = find(qual_' num2str(IYear) Filler num2str(IStart) 'd1.qual < IQual);'])  % Night
    eval(['SSTd1 = sst_' num2str(IYear) Filler num2str(IStart) 'd1.sst;'])

    eval(['nnd3 = find(qual_' num2str(IYear) Filler num2str(IStart) 'd3.qual < IQual);'])  % Day
    eval(['SSTd3 = sst_' num2str(IYear) Filler num2str(IStart) 'd3.sst;'])

    SSTd1(nnd1) = nan;
    SSTd3(nnd3) = nan;
else
    eval(['SSTd1 = sst_' num2str(IYear) Filler num2str(IStart) 'd1.sst;'])  % Night
    eval(['SSTd3 = sst_' num2str(IYear) Filler num2str(IStart) 'd3.sst;'])  % Day
end

for i=IStart:IEnd
    if IStart >= 100
        Filleri = '';
    else if IStart >= 10
            Filleri = '0';
        else
            Filleri = '00';
        end
    end

    if IQual > 0
        eval(['nn = find(qual_' num2str(IYear) Filleri num2str(i) 'd1.qual < IQual);'])  % Night
        eval(['SSTTemp1 = sst_' num2str(IYear) Filleri num2str(i) 'd1.sst;'])
        SSTTemp1(nn) = nan;

        eval(['nn = find(qual_' num2str(IYear) Filleri num2str(i) 'd3.qual < IQual);'])  % Day
        eval(['SSTTemp3 = sst_' num2str(IYear) Filleri num2str(i) 'd3.sst;'])
        SSTTemp3(nn) = nan;
    else
        eval(['SSTTemp1 = sst_' num2str(IYear) Filleri num2str(i) 'd1.sst;'])
        eval(['SSTTemp3 = sst_' num2str(IYear) Filleri num2str(i) 'd3.sst;'])
    end

    SSTd1 = max( SSTd1, SSTTemp1);
    SSTd3 = max( SSTd3, SSTTemp3);

end

% Plot night
figure(ImagePlane)
eval(['image( sst_' num2str(IYear) Filler num2str(IStart) 'd1.lon, sst_' num2str(IYear) Filler num2str(IStart) 'd1.lat, SSTd1*64/' num2str(MaxTemperature) ')'])
set(gca,'ydir','normal')
title('SST for Night Passes Composited')

%Plot day
figure(ImagePlane+1)
eval(['image( sst_' num2str(IYear) Filler num2str(IStart) 'd3.lon, sst_' num2str(IYear) Filler num2str(IStart) 'd3.lat, SSTd3*64/' num2str(MaxTemperature) ')'])
set(gca,'ydir','normal')
title('SST for Day Passes Composited')

%Plot day and night together
SST= max(SSTd1,SSTd3);
figure(ImagePlane+2)
eval(['image( sst_' num2str(IYear) Filler num2str(IStart) 'd3.lon, sst_' num2str(IYear) Filler num2str(IStart) 'd3.lat, SST*64/' num2str(MaxTemperature) ')'])
set(gca,'ydir','normal')
title('SST for Day and Night Passes Composited')
