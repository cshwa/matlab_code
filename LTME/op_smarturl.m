function URL=op_smarturl(URL,FIELD,kT)
if isstruct(URL) 
    %each field has its own URL
    URL=eval(['URL.' FIELD]);
end
if iscell(URL)
    % URL is different for each time frame
    URL=URL{kT};
end
%if strfind(URL,'URL') & strfind(URL,'=') % if string is a metastring URL='URL=something' 
%    eval(URL); 
%end
