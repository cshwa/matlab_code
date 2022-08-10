nx = 51;
ny = 51;
ntot = 1000;
nout = 5;

fid = fopen('eta.dat','r');
height = fread(fid,nx*ny*(ntot/nout),'single');
eta=reshape(height,[nx ny ntot/nout]);
fclose(fid);

fid = fopen('eta_adj.dat','r');
height = fread(fid,nx*ny*(ntot/nout),'single');
eta_adj=reshape(height,[nx ny ntot/nout]);
fclose(fid);

set(gcf,'position',[100 100 1000 800]);

    vidObj = VideoWriter('eta.avi');
    open(vidObj);
 
    % Create an animation.
    axis tight
    set(gca,'nextplot','replacechildren');

    [x,y]=meshgrid([1:1:51]);
    
    for k = 1:ntot/nout
       clf;
       contour3(x,y,eta_adj(:,:,k));
       surface(x,y,eta_adj(:,:,k),'EdgeColor',[0.8,0.8,0.8],'FaceColor','none');
       axis([1 50 1 50 -0.1 0.1]);
       view([40 40 20]);
%        contour3(eta(:,:,k),'fill','on');
%        caxis([-0.5 0.5]);
%        colorbar;
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
    end
  
    % Close the file.
    close(vidObj);