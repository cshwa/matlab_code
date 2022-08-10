function arrow_color(hq1,headsize)

hkid = get(hq1,'children');
XX = get(hkid(1),'XData');
YY = get(hkid(1),'YData');
cmap = jet(361); %colormap

for ii = 1:3:length(XX)-1
    headWidth = headsize * sqrt((XX(ii+1)-XX(ii)).^2 + (YY(ii+1)-YY(ii)).^2); % set the headWidth, function of length of arrow
    angled = floor(atan2(YY(ii+1)-YY(ii),XX(ii+1)-XX(ii))*180/pi) + 181; %get the angle
    if isnan(angled)==0
    ah = annotation('arrow',...
        'Color', cmap(angled,:),...
        'headStyle','cback1','HeadLength',30,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[XX(ii) YY(ii) XX(ii+1)-XX(ii) YY(ii+1)-YY(ii)]);
    end

end

