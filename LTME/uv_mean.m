clearvars; close all;
mags=[10,8,7];
% angs=[95,85,75];
angs=[0,180,75];
u = mags.*cosd(angs);
v = mags.*sind(angs);

%% component-wise mean
u_comp_m = mean(u);
v_comp_m = mean(v);
%% trigonometric mean
u_trig_m = mean(mags)*mean(cosd(angs));
v_trig_m = mean(mags)*mean(sind(angs));
%% circular mean
mag_m = mean(mags);
ang_m = circular_mean(angs);
u_circ_m = mag_m*cosd(ang_m);
v_circ_m = mag_m*sind(ang_m);
%% figure
figure; hold on; grid on; box on;
scatter(u,v,50,'b'); leg_str = {'original'};
scatter(u_comp_m,v_comp_m,50,'r'); leg_str = [leg_str(:)',{'component-wise mean'}];
scatter(u_trig_m,v_trig_m,50,'g'); leg_str = [leg_str(:)',{'trigonometric mean'}];
scatter(u_circ_m,v_circ_m,50,'k'); leg_str = [leg_str(:)',{'circular mean'}];
axis equal; xlim([-10,10]); ylim([-10,10]);
legend(leg_str,'location','eastoutside');

%% function
function thbar = circular_mean( angs )
    sbar = mean( sind(angs), 'all' );
    cbar = mean( cosd(angs), 'all' );
    thbar = atand( sbar/cbar );
    if cbar < 0
        thbar = thbar + 180;
    elseif sbar < 0
        thbar = thbar + 360;
    end
end