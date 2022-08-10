function varargout = op_timemessage(OrbitType, MoreThanDaily)
% Generate time explanation for satellite-derived data.
%
% op_template
% h = op_template(OrbitType,MoreThanDaily)
%
% INPUT(S):
%
%   OrbitType -- 'polar','goes','pathfinder1km','seawinds'
%   MoreThanDaily -- ???
%
%   if MoreThanDaily == 1 the following text is added:
%
%   'Averaging periods other than daily are generated by averaging the 
%    daily ascending (or descending) fields over the indicated interval. 
%    For example, if ascending fields and 8 day Temporal Averages were 
%    selected, the returned fields would consist of one field every 8 days 
%    obtained by averaging 8 consecutive daily ascending fields.'
%
% OUTPUT(S):
%
%   h -- (optional) handle to message box
%
% EXAMPLE:
%
% OPeNDAP Science Team
% Copyright 2008
% $Version 2.0.2$

%==========================================================================
% Peter Cornillon
%
% REVISION HISTORY:
% 2008/02/15 1.0.0 (beta) created, pcornillon
% 2008/02/28 1.0.1 (beta) edited, ceb
% 2008/03/01 2.0.0 released
% 2008/03/03 2.0.1 fixed typo on 'seawinds' section, ceb
% 2008/03/18 2.0.2 edited "help" section, ceb
%==========================================================================

CreateMode.Interpreter = 'Tex';
CreateMode.WindowStyle = 'Modal';

OrbitType = lower(OrbitType);

switch OrbitType

    case 'polar'

        TimeMessageTxt = cell(7,1);
        if MoreThanDaily
            TimeMessageTxt = cell(9,1);
        end
        ii = 1;
        TimeMessageTxt{ii} = ' '; ii = ii+1;
        TimeMessageTxt{ii} = ['Polar orbiting satellites circle the Earth approximately', ...
            ' 14 times per day. For one half of each orbit, the satellite is moving from', ...
            ' south to north (ascending) and for the other half from north to south', ...
            ' (descending). For sun-synchronous satellites the ascending portions of', ...
            ' the orbits occur at the same local sun time (LST) everywhere on Earth, e.g.,', ...
            ' a 2PM (ascending) pass would fly overhead 2 hours after the sun had reached its', ...
            ' highest point in the sky.', ...
            ' Descending portions of the orbits occur approximately 12 hours later.\newline']; ii = ii+1;
        TimeMessageTxt{ii} = ' '; ii = ii+1;
        TimeMessageTxt{ii} = ['Global data products are generally developed by combining all', ...
            ' ascending portions of the passes from a day together to form one map and', ...
            ' all descending portions of the passes to form a second map. Each of these is',...
            ' referred to as a daily average because it consists of data collected over', ...
            ' an approximately 24 hour period. The time of the field is associated with', ...
            ' the local sun time; e.g., the 2PM field, also referred to as the ascending',...
            ' orbit or pass. The other field would be referred', ...
            ' to as the 2AM field or the descending orbit or pass. \newline'];ii = ii+1;
        TimeMessageTxt{ii} = ' '; ii = ii+1;
        TimeMessageTxt{ii} = ['BOTTOM LINE: Global fields are', ...
            ' daily averages separated by approximately 12 hours, alternating', ...
            ' between ascending and descending fields.\newline']; ii = ii+1;
        TimeMessageTxt{ii} = ' '; ii = ii+1;

        if MoreThanDaily
            TimeMessageTxt{ii} = ['Averaging periods other than daily are generated by', ...
                ' averaging the daily ascending (or descending) fields over the', ...
                ' indicated interval. For example, if ascending fields and 8 day', ...
                ' ''Temporal Averages'' were selected, the returned fields would consist of one', ...
                ' field every 8 days obtained by averaging 8 consecutive daily', ...
                ' ascending fields.']; ii = ii+1;
            TimeMessageTxt{ii} = ' '; ii = ii+1;
        end

        h=msgbox( TimeMessageTxt,'Temporal Sampling and Averaging',CreateMode);

    case 'goes'

        TimeMessageTxt = cell(5,1);
        ii = 1;
        TimeMessageTxt{ii} = ' '; ii = ii+1;
        TimeMessageTxt{ii} = ['The orbits of geosynchronous satellites lie in the', ...
            ' Earth''s equatorial plane and the height of the satellite above the Earth is', ...
            ' chosen so that these satellites circle the Earth once per day;', ...
            ' i.e., they remain above the same geographic location on the', ...
            ' equator through time. As a result, geosynchronous satellites', ...
            ' see only a portion', ...
            ' of the Earth''s surface and the spatial resolution descreases away', ...
            ' from the satellite subpoint. These satellites make up for their', ...
            ' coarse spatial resolution with high temporal resolution.', ...
            ' They scan the portion of the Earth''s', ...
            ' surface that they see every 30 minutes. \newline']; ii = ii+1;
        TimeMessageTxt{ii} = ' '; ii = ii+1;
        TimeMessageTxt{ii} = ['For this data set an average SST field is generated for',...
            ' each day in the specified temporal interval. Averaging is performed', ...
            ' over all fields falling within a temporal range specified by',...
            ' the selected ''Temporal Averages''. A daily ''Temporal Averages'''....
            ' returns one field per day averaged over the 48 fields of that day.',...
            ' An 8 day ''Temporal Averages'' returns one field per day averaged',...
            ' over the 388 fields in the 8 day interval centered on that day. \newline']; ii = ii+1;
        %TimeMessageTxt{ii} = ' '; ii = ii+1;

        h = msgbox(TimeMessageTxt,'GOES Temporal Sampling and Averaging',CreateMode);
        

	case 'pathfinder1km'

        TimeMessageTxt = cell(5,1);

        TimeMessageTxt{1} = ' ';
        TimeMessageTxt{2} = ['Polar orbiting satellites circle the Earth approximately', ...
            ' 14 times per day. For one half of each orbit, the satellite is moving from', ...
            ' south to north (ascending) and for the other half from north to south', ...
            ' (descending). For sun-synchronous satellites the ascending portions of', ...
            ' the orbits occur at the same local sun time (LST) everywhere on Earth, e.g.,', ...
            ' a 2PM (ascending) pass would fly overhead 2 hours after the sun had reached its', ...
            ' highest point in the sky.', ...
            ' Descending portions of the orbits occur approximately 12 hours later.\newline'];
        TimeMessageTxt{3} = ' ';
        TimeMessageTxt{4} = ['For these data sets, only small portions of each orbit,',...
            ' corresponding to the data collected while the satellite was in view of',...
            ' a ground receiving station are available. These segments are approximately',...
            ' 11 minutes or 4000km long. All of the orbits with substantial coverage',...
            ' of ocean waters were collected. The actual number of orbits per day',...
            ' varies between 3 and 5 depending on the location of the orbits for that day.',...
            ' Consecutive ascending (or descending) segments are separated by 90 minutes.',...
            ' On average there will be 2 ascending and 2 descending passes, the 2',...
            ' ascending passes separated by 90 minutes, the 2nd ascending and 1st descending',...
            ' passes separated by about 11 hours and the 2 descending passes separated',...
            ' by 90 minutes. \newline'];
        %TimeMessageTxt{5} = ' ';
        
        h=msgbox( TimeMessageTxt,'Pathfinder 1km Temporal Sampling',CreateMode);
        
    case 'seawinds'

        TimeMessageTxt = cell(9,1);

        TimeMessageTxt{1} = ' ';
        TimeMessageTxt{2} = ['Polar orbiting satellites circle the Earth approximately', ...
            ' 14 times per day. For one half of each orbit, the satellite is moving from', ...
            ' south to north (ascending) and for the other half from north to south', ...
            ' (descending). For sun-synchronous satellites the ascending portions of', ...
            ' the orbits occur at the same local sun time (LST) everywhere on Earth, e.g.,', ...
            ' a 6PM (descending) pass would fly overhead 6 hours after the sun had reached its', ...
            ' highest point in the sky.', ...
            ' Descending portions of the orbits occur approximately 12 hours later.\newline'];
        TimeMessageTxt{3} = ' ';
        TimeMessageTxt{4} = ['Global data products are generally developed by combining all', ...
            ' ascending portions of the passes from a day together to form one map and', ...
            ' all descending portions of the passes to form a second map. Each of these is',...
            ' referred to as a daily average because it consists of data collected over', ...
            ' an approximately 24 hour period. The time of the field is associated with', ...
            ' the local sun time; e.g., the 6PM field, aslo referred to as the descending',...
            ' orbit or pass. The other field would be referred', ...
            ' to as the 6AM field or the ascending orbit or pass. \newline'];
        TimeMessageTxt{5} = ' ';
        TimeMessageTxt{6} = ['BOTTOM LINE: Global fields are', ...
            ' daily averages separated by approximately 12 hours, alternating', ...
            ' between ascending and descending fields.\newline'];
        TimeMessageTxt{7} = ' ';
        TimeMessageTxt{8} = ['Examples: Selecting one pass time, say 6:00 a.m., every 3 time steps, will return',...
            ' one ascending field for every third day in the specified temporal interval. Selecting',...
            ' both passes and every time step will return all of the fields in the',...
            ' specified interval. Finally selecting both passes every 3 time steps will return',...
            ' two fields for every third day, one 6 a.m. and one 6 p.m. These fields',...
            ' are not averaged over the 3 day interval, they are sampled every 3 days. .\newline'];
        TimeMessageTxt{9} = ' ';

        h=msgbox( TimeMessageTxt,'Temporal Sampling and Averaging',CreateMode);
end

if nargout
    varargout{1} = h;
end

return
