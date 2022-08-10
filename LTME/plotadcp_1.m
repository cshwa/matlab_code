% plot the adcp data as direction arrows.  Loads the .mat files created by
% adcp.m.  Must be run from the same directory the files are in.
%
% Bill Scuba   Oct. 15, 2003
% Scripps Institution of Oceanography, San Diego, CA, USA
% wscuba@ucsd.edu

clear dirstruct

figure
hold on
dirstruct=dir('*.mat')
for yi=1:length(dirstruct)
    clear ping 
    load(dirstruct(yi).name)
	for xi=1:50:length(ping) % plots every ~50th ping
		quiver(ping(xi).lon,ping(xi).lat,ping(xi).adcpdata(1,4),ping(xi).adcpdata(1,5),.002)
	end
    fprintf('u = %g \n',sqrt(ping(xi).adcpdata(1,4)^2 +ping(xi).adcpdata(1,5)^2))
end
hold off
