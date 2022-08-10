function op_plotimage(FigNo,StructureName,NoCoast)
% op_plotimage(FigNo,StructureName,NoCoast)
%
% Plots images of data obtained with a downloaded from one of the OPeNDAP 
% GUIs. The user is prompted for input based on the data to be plotted. The
% function will plot images of scalar arrays or vectors of vector arrays.
% Vector plotting is performed under user control. Vectors may be plotted
% over scalar fields, again under user control after submitting the plot
% request. The fields need not be obtained from the same GUI. If the
% coordinate system is rectangular in longitude and latitude, the scalar
% plots are performed with imagesc. If the coordinate system is not
% rectangular in longitude and latitude, scalar plots are performed with
% pcolor.
%
% INPUT(S):
%
%   FigNo -- number of the figure in which image will be plotted.
%   StructureName -- name of structure in user's workspace.
%   NoCoast -- (optional) If present will not try to plot coastline.
%
% OUTPUT(S):
%
%   [none] --
%
% NOTE: If the field has more than 2 dimensions, the user is prompted for the
% dimensions to be suppresed.
%
% EXAMPLE: Acquire GOES data and SeaWinds data either via the respective
% GUIs or as follows and plot the vectors on the SST data.
%
% 1) Generate a structure called request_pathfinder1km with the following content:
% 
% request_pathfinder = 
% 
%            Fields: {'sst_dec'  'sst_raw'}
%       Coordinates: {'latitude'  'longitude'  'time'}
%            Passes: {'night'  'day'}
%             DATE1: '19940401'
%             DATE2: '19940401'
%              LAT1: 0
%              LAT2: 66.00
%              LON1: -164.00
%              LON2: 1.00
%          DATEINCR: 1.00
%           LATINCR: 8.00
%           LONINCR: 8.00
%           DIRNAME: ''
%          FNPREFIX: ''
%     SaveWorkspace: 'y'
%         SaveFiles: 'n'
%          SaveMode: 'opendap'
%     DataSetBranch: 'Pathfinder1km_NWA'
%       RequestDate: '2009-09-15 13:32:29'
%         iTimeIncr: 1.00
%         Variables: {'sst_dec'  'sst_raw'  'latitude'  'longitude'  'time'}
%
%  2) Get Pathfinder data:
%
%    >> get_pathfinder1km(request_pathfinder)
%
%  This will generate three structures: opendap_0001 - opendap_0003 with
%  the pathfinder data of interest.
%
%  3) Generate a structure called request_hycom:
%
% request_hycom = 
% 
%            Fields: {'mixed_layer_u_velocity'  'mixed_layer_v_velocity'}
%       Coordinates: {'latitude'  'longitude'  'time'}
%             DATE1: '20040401'
%             DATE2: '20040401'
%              LAT1: 30.00
%              LAT2: 45.00
%              LON1: -80.00
%              LON2: -50.00
%            DEPTH1: []
%            DEPTH2: []
%          DATEINCR: 1.00
%           LATINCR: 1.00
%           LONINCR: 1.00
%         LEVELINCR: 1.00
%           DIRNAME: ''
%          FNPREFIX: ''
%     SaveWorkspace: 'y'
%         SaveFiles: 'n'
%          SaveMode: 'opendap'
%     DataSetBranch: 'HYCOM_Global_Hindcast'
%       RequestDate: '2009-09-15 15:21:20'
%         iTimeIncr: 1.00
%         Variables: {'mixed_layer_u_velocity'  'mixed_layer_v_velocity'  'latitude'  'longitude'  'time'}
%
%  This will generate 1 structures with the Hycom data.
%
% > Skip until we get the merging function working again.
% > 4) merge the opendap_nnn fields using the op_mergingtool command. This
% > command will merge opendap_0001-0003 into one structure called goes and
% > opendap_0004-0009 into a second structure called qscat. You could have
% > plotted directory from the opendap structures, but for this example, we
% > will plot from the merged structures.
%
%  5) Plot the Pathfinder SST data:
%
%   >> op_plotimage(1,opendap_0002)
%
%  sst is a (latitude,longitude,time) array.
%  op_plotimage will only plot (longitude,latitude) arrays so you must select a single value in time.
%  Acceptable values of time are:
%  1: 2008-03-05 12:00:00
%  2: 2008-03-06 12:00:00
%  3: 2008-03-07 12:00:00
%  Choose one of the above values [1-3]: 2
%
% This will plot the GOES data in figure 1.
%
%  6) Plot the SeaWinds data:
%
%  op_plotimage(1,qscat)
%  Looks like wind_u and wind_v are vector components, would you like to plot them with quiver [y/n]: y
%  Do you want to plot vectors on the magnitude of the vector? [y/n]: n
%  Enter scale to use when vectors are plotted - default is 1.0? 3
%  Enter the color in which the vectors are to be plotted [w,r,g,b,k,...] default is w? k
%  wind_uwind_v is a (latitude,longitude,time) array.
%  op_plotimage will only plot (longitude,latitude) arrays so you must select a single value in time.
%  Acceptable values of time are:
%  1: 2008-03-05 06:00:00
%  2: 2008-03-05 18:00:00
%  3: 2008-03-06 06:00:00
%  4: 2008-03-06 18:00:00
%  5: 2008-03-07 06:00:00
%  6: 2008-03-07 18:00:00
%  Choose one of the above values [1-6]: 4
%  Figure number 1 already exists. Would you like to overplot [y/n] Default n? y
%
%   This will plot the SeaWinds vectors in black on the GOES SST image.
%
% OPeNDAP Science Team
% Copyright 2007, 2008
% $Version 2.0.1$

%==========================================================================
% Peter Cornillon
%
% REVISION HISTORY:
% 2006/12/01 0.x.x created
% 2007/12/01 0.x.x modified for new data structure
% 2008/02/28 2.0.0 prepare for distribution
% 2008/03/18 2.0.1 updated help, fixed coastal plotting
% 2008/03/18 2.1.0 Added quiver capability, rewrote scalar plotting, fixed
% to handle pcolor, added help and example.
%==========================================================================

% Cell array of metadata field names; fields that will not be plotted.

MetadataNames{1} = 'latitude';
MetadataNames{2} = 'longitude';
MetadataNames{3} = 'time';
MetadataNames{4} = 'user_friendly_time';
MetadataNames{5} = 'dimensions';
MetadataNames{6} = 'variables';
MetadataNames{7} = 'attributes';

NumArguments = nargin;
if NumArguments < 2
    disp(' ')
    disp(' Must supply either 2 or 3 arguments for this function.')
    disp('    First argment is the figure number in which to plot the results.')
    disp('    Second argument is the name of the structure containing the data to be plotted.')
    disp(' Try it again.')
    disp(' ')
    return
end

% Get field names of the structure.

FieldNames = fieldnames( StructureName);   % Get field names first.
FieldNameslc = lower(FieldNames);   % Shift to lower case to avoid case problems in the comparison.
NumFieldNames = length(FieldNameslc);


%% Now get list of plottable data; actually exclude non-plottable fields from list.

iSave = [];
for iNames=1:NumFieldNames
    isave = 1;
    for jNames=1:length(MetadataNames)
        if ~isempty(findstr(FieldNameslc{iNames},MetadataNames{jNames}))
            isave = 0;
        end
    end

    % This variable should be plottable; add it to the list.

    if isave == 1
        iSave =[iSave iNames];
    end
end

NumArguments = length(iSave);


% Make sure that a variable to plot has been selected.

QuiverPlot = 'n';  % Need to set this. May reset it below.

if NumArguments == 1
    VariableToPlot = iSave;
else

    % Here if more than one plottable variable. First check to see if any of the pairs are u,v

    clear iSaveU iSaveV
    for i=iSave
        if strfind(FieldNames{i},'_u')
            iSaveU = i;
        end
        if strfind(FieldNames{i},'_v')
            iSaveV = i;
        end
    end

    % Figure out if the variable list contains u,v vectors.
    % clear iSaveU iSaveV; QuiverPlot = 'n';   % Temporary to skip quiver plotting stuff.

    QuiverPlot = 'n';  % Need a value of QuiverPlot for later
    Vector = 0;
    if exist('iSaveU') && exist('iSaveV')
        Vector = 1;
        QuiverPlot = input(['Looks like ' FieldNames{iSaveU} ' and ' FieldNames{iSaveV} ...
            ' are vector components, would you like to plot them with quiver [y/n]: '], 's');
        if isempty(QuiverPlot) % Default is yes
            QuiverPlot = 'y';
        end
        if length(QuiverPlot) == 2  % If QuiverPlot is 2 characters assume user entered no and change to n
            QuiverPlot = 'n';
        end
        if length(QuiverPlot) == 3  % If QuiverPlot is 3 characters assume user entered yes and change to y
            QuiverPlot = 'y';
        end

        if QuiverPlot == 'y'   
            PlotSpeed = input('Do you want to plot vectors on the magnitude of the vector? [y/n]: ', 's');
            VectorLength = input('Enter scale to use when vectors are plotted - default is 1.0? ');
            if isempty(VectorLength)
                VectorLength = 1.0;
            end
            VectorColor = input('Enter the color in which the vectors are to be plotted [w,r,g,b,k,...] default is w? ', 's');
            if isempty(VectorColor)
                VectorColor = 'w';
            end
        end
    end

    if QuiverPlot == 'y'
        VariableToPlot = [iSaveU iSaveV];
    else
        if NumArguments == 2 && Vector == 1
            VariableToPlot = input(['OK which of these variables would you like to plot [1 for ' ...
                FieldNames{min(iSaveU,iSaveV)} ' or 2 for ' FieldNames{max(iSaveU,iSaveV)} ']: ']);
        else
            disp(' Plottable variables available for this array are')
            for i=iSave
                disp(['   ' num2str(i) ': ' FieldNames{i}])
            end
            VariableToPlot = input(['Choose one of these [1-' iSave ']: '])
        end
    end
end

% Get the name of the variable to plot - for now will not deal with vectors.

eval(['PlotVariable = cell2mat(FieldNames(' num2str(VariableToPlot(1)) '));'])    % Get the name of the element to plot

if QuiverPlot == 'y'
    eval(['PlotVariable2 = cell2mat(FieldNames(' num2str(VariableToPlot(2)) '));'])    % Get the name of the element to plot
end

% eval([ 'StructureName.' PlotVariable ' = squeeze(StructureName.' PlotVariable ');'])

%% Get the name used for the dimension(s) of latitude and longitude

LatIs1D = 1;
if length(StructureName.Variables.latitude) == 2
    LatIs1D = 0;
end

LatDim = StructureName.Variables.latitude{1};
LonDim  = StructureName.Variables.longitude{1};
if LatIs1D == 0
%    LatDim(2) = StructureName.Variables.latitude{2};
    LonDim  = StructureName.Variables.longitude{2};
end

%% Get the dimension names of the variable to plot.

eval(['Dimensionality = size(StructureName.' PlotVariable ');'])

NDimensions = length(Dimensionality);

for iDimension=1:NDimensions   % Get names of dimensions of variable to plot.
%    eval(['DimensionNames{' num2str(iDimension) '} = StructureName.Dimensions.' PlotVariable '_variables{' num2str(iDimension) '};'])
    eval(['DimensionNames{' num2str(iDimension) '} = StructureName.Variables.' PlotVariable '{' num2str(iDimension) '};'])
end

for iDimension=1:NDimensions
    if iDimension == 1
        VarDimensions = ['(' DimensionNames{iDimension}];
    else
        VarDimensions = [VarDimensions ',' DimensionNames{iDimension}];
    end
end
VarDimensions = [VarDimensions ')'];

% And associate them with lat and lon dimensions.

iLat = 0; 
iLon = 0;

for iDimension=1:NDimensions
    if findstr(DimensionNames{iDimension},LatDim)
        iLat = iDimension;
    end
    if findstr(DimensionNames{iDimension},LonDim)
        iLon = iDimension;
    end
end

% Now associate these dimension names with actual variables.

for iNames=1:NumFieldNames
%    if strcmp(FieldNameslc{iNames},DimensionNames{iLat})
    if strcmp(FieldNameslc{iNames},'latitude')
        iLatitude = iNames;
    end
%    if strcmp(FieldNameslc{iNames},DimensionNames{iLon})
    if strcmp(FieldNameslc{iNames},'longitude')
        iLongitude = iNames;
    end
end

%% This function only plots longitude, latitude arrays, abort if not present.


if iLat * iLon == 0
    disp('op_plotimage will only plot arrays in longitude and latitude.')
    disp(['Dimensions of ' PlotVariable ' are: ' VarDimensions])
    return
end

%% If more than two dimensions let the user know that they will have to subset and on which dimensions,

% Construct string to print out for non-lat,lon dimensions if more than 3 dimensions
% Also get the list of non-lat,lon dimensions ==> residual dimensions

if NDimensions >= 3

    iRD = 0;

    for iDimension=1:NDimensions

        if ~strcmp(DimensionNames{iDimension},'Y') && ~strcmp(DimensionNames{iDimension},'X')

            iRD = iRD + 1;

            if iRD == 1
                ResidualDimensions =  DimensionNames{iDimension};
            elseif iRD < NDimensions-2
                ResidualDimensions = [ResidualDimensions ', ' DimensionNames{iDimension}];
            else
                ResidualDimensions = [ResidualDimensions ' and ' DimensionNames{iDimension}];
            end

        end
    end

    disp([FieldNames{VariableToPlot} ' is a ' VarDimensions ' array.'])
    disp(['op_plotimage will only plot (longitude,latitude) arrays so you must select a single value in ' ResidualDimensions '.'])

end

%% Now we get the slice of the residual dimensions and plot

Constraint = ['('];
iFirst = 1;
for iDimension=1:NDimensions
%   if strcmp(DimensionNames{iDimension},'latitude') | strcmp(DimensionNames{iDimension},'longitude')
    if iLat == iDimension | iLon == iDimension
        if iFirst == 1
            Constraint = [Constraint ':'];
            iFirst = 0;
        else
            Constraint = [Constraint ',:'];
        end
    else

        % Now get the slice to plot for this dimension

        if Dimensionality(iDimension) < 10  % If less than 10 values list them. Otherwise just ask for number.
            disp(['Acceptable values of ' DimensionNames{iDimension} ' are:'])
            for iValues=1:Dimensionality(iDimension)
                if strcmp(DimensionNames{iDimension},'time')
                    disp([num2str(iValues) ': ' StructureName.metadata.user_friendly_time{iValues}])
                else
                    eval(['xx = StructureName.' DimensionNames{iDimension} '(iValues);'])
                    disp([num2str(iValues) ': ' num2str(xx)])
                end
            end
            ElementToPlot = input(['Choose one of the above values [1-' num2str(Dimensionality(iDimension)) ']: ']);
        else
            ElementToPlot =   input(['There are ' num2str(Dimensionality(iDimension)) ' of ' DimensionNames{iDimension} ...
                '. Specify which value of this dimension you would like to plot: ']);
        end

        if iFirst == 1
            Constraint = strcat(Constraint,num2str(ElementToPlot));
            iFirst = 0;
        else
            Constraint = strcat(Constraint,',',num2str(ElementToPlot));
        end
    end
end

Constraint = strcat(Constraint,')');

%% Now plot; check to see if plot already exists.

% If the plot already exists and if this is a quiver plot ask if they user
% wants to overlay the plots, an image plot on top of another one doesn't
% make a lot of sense. 

OverPlot = 'n';
if ishandle(FigNo) && QuiverPlot == 'y'&& PlotSpeed == 'n'
    OverPlot = input(['Figure number ' num2str(FigNo) ' already exists. Would you like to overplot [y/n] Default n? '], 's');
    if isempty(OverPlot)
        OverPlot = 'n';
    end
end

if OverPlot == 'n'
    eval(['figure(' num2str(FigNo) ')'])
    hold off
else
    eval(['figure(' num2str(FigNo) ')'])
    hold on
end

%% Is this a regular longitude, latitude grid, or do I have to use pcolor?

% I will have to use pcolor in one of two cases: the lat and or lon
% vectors are not linear or they are two dimensional. Check to see of 2d.

LinearFlag = 0;
if LatIs1D == 1
    LatStep = diff(StructureName.latitude);
    LonStep = diff(StructureName.longitude);
    LinearFlag = (max(LatStep)-min(LatStep)) < .000001 && (max(LonStep)-min(LonStep)) < .000001;
end

if LinearFlag == 1 && LatIs1D == 1

    ContinentalOutlineColor = 'w';
   
    if QuiverPlot == 'n'
        eval(['imagesc( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
            ', squeeze(StructureName.' PlotVariable Constraint '))'])
        colorbar
        set(gca,'ydir','normal')
    else
        if PlotSpeed == 'y'
            eval(['u = squeeze(StructureName. ' PlotVariable Constraint ');'])
            eval(['v = squeeze(StructureName. ' PlotVariable2 Constraint ');'])
            Speed = sqrt( u .* u + v.* v);
            eval(['imagesc( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
                ', Speed)'])
            colorbar
            set(gca,'ydir','normal')
            hold on
        end
        
        eval(['quiver( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
            ', squeeze(StructureName.' PlotVariable Constraint ...
            '), squeeze(StructureName.' PlotVariable2 Constraint '),' num2str(VectorLength) ', ''' VectorColor ''')'])
 
        % Need to set axes if no image and coastline is to be plotted.
        
        LonMin = min(StructureName.longitude);
        LonMax = max(StructureName.longitude);
        LatMin = min(StructureName.latitude);
        LatMax = max(StructureName.latitude);
        axis([LonMin LonMax LatMin LatMax])
        
        % Set color for continental outline.
        
        ContinentalOutlineColor = 'b';

    end

else % Have to use pcolor.

    ContinentalOutlineColor = 'k';

    if QuiverPlot == 'n'
        eval(['pcolor( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
            ', squeeze(StructureName.' PlotVariable Constraint '))'])
        colorbar
        shading flat
        set(gca,'ydir','normal')
    else
        if PlotSpeed == 'y'
            eval(['u = squeeze(StructureName. ' PlotVariable Constraint ');'])
            eval(['v = squeeze(StructureName. ' PlotVariable2 Constraint ');'])
            Speed = sqrt( u .* u + v.* v);
            eval(['pcolor( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
                ', Speed)'])
            shading flat
            colorbar
            set(gca,'ydir','normal')
            hold on
            eval(['quiver( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
                ', squeeze(StructureName.' PlotVariable Constraint ...
                '), squeeze(StructureName.' PlotVariable2 Constraint '),' num2str(VectorLength) ', ''' VectorColor ''')'])

        else
            
            % First need to determine if the longitude values of the
            % vectors fall in the range of the plot. If not will shifting
            % them by 360 accomplish this? If so shift them to plot.
            
            LonRange = get(gca,'xlim');
            LonRangeQuiver(1) = min(StructureName.longitude(:));
            LonRangeQuiver(2) = max(StructureName.longitude(:));
            
            LatRange = get(gca,'ylim');
            
            Longitude = StructureName.longitude;
            if min(LonRange) < 0 && max(LonRangeQuiver) > 180
                nn = find(StructureName.longitude > 180);
                Longitude(nn) = StructureName.longitude(nn) - 360;
            end
            
            if max(LonRange) > 180 & min(LonRangeQuiver) < 0
                nn = find(StructureName.longitude < 0); 
                Longitude(nn) = StructureName.longitude(nn) + 360;
            end
            
%            eval(['quiver( StructureName.' FieldNames{iLongitude} ', StructureName.' FieldNames{iLatitude} ...
            eval(['quiver( Longitude, StructureName.' FieldNames{iLatitude} ...
                ', squeeze(StructureName.' PlotVariable Constraint ...
                '), squeeze(StructureName.' PlotVariable2 Constraint '),' num2str(VectorLength) ', ''' VectorColor ''')'])
%            LonMin = min(StructureName.longitude);
%            LonMax = max(StructureName.longitude);
%            LatMin = min(StructureName.latitude);
%            LatMax = max(StructureName.latitude);
            LonMin = min(min(Longitude(:)), LonRange(1));
            LonMax = max(max(Longitude(:)), LonRange(2));
            LatMin = min(min(StructureName.latitude(:)), LatRange(1));
            LatMax = max(max(StructureName.latitude(:)), LatRange(2));
            
            axis([LonMin LonMax LatMin LatMax])

            % Set color for continental outline.

            ContinentalOutlineColor = 'b';
            
        end

    end
end


%% Annotate the plot, add coastline if available.


% Now plot the coast. I think that coast is in image toolbox, so if this
% toolbox is not present, disable the following.

ToolBoxes = ver;
for i=1:length(ToolBoxes)
    if strcmp(ToolBoxes(i).Name,'Mapping Toolbox') == 1   % Mapping toolbox present?

        load coast   % Load the coast.

        eval(['NeedToShift = max(StructureName.' FieldNames{iLongitude} ') > 180; '])

        if NeedToShift   % This section shifts coastline from -180 to 180 to 0 to 360.
            [Lon Lat] = op_shiftcoast(long,lat);
        else
            Lon = long;
            Lat = lat;
        end
        if ~exist('NoCoast')
            hold on
            plot(Lon,Lat,ContinentalOutlineColor)
        end
    end
end


return
