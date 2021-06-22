%% PLOT IPP video
% Author :Fabio Manta 
% Last update: Dec-2020

      
% load coastsSEA.dat
% llon=coastsSEA(:,1);
% llat=coastsSEA(:,2);

id_staz = (1:length(staz));
SatNum = [18 25 26 29 31];
filename = 'movie_2.gif';
initial_time = 15; % in hours
dt= 30;
timeStep = (0:dt:3600*2);
delta_deg =10;

figure;

for tim = 1:size(timeStep,2) %time loop
%     plot(llon,llat,'k')
        
    hold on
    
    for jj = 1:size(id_staz,2)
         
        for ii = 1:length(SatNum)  
%     hold on    
         % plot all station tecu
            tmpTIME = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,1)',SatNum(ii)))); 
            tmpLAT = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,2)',SatNum(ii))));  
            tmpLON = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,3)',SatNum(ii))));
            
            tmpTecu = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,8)',SatNum(ii))));                       
    
            startingPosition = find(tmpTIME <= initial_time);
%             startingPosition = find(tmpTIME <= EventTimeHours+0);
            
                    if ~isempty (startingPosition) 
                        startingPosition = startingPosition(end);
                        
                        if startingPosition+tim-1<=(length(tmpLON))
             
                        cplot(tmpLON(startingPosition+(tim)-1),tmpLAT(startingPosition+(tim)-1)...
                        ,tmpTecu(startingPosition+(tim)-1),'o','markerfacecolor','flat','MarkerSize',5);
                        
                        end
                     
                        
                    end

        end
    end
    
    [Y,M,D,H,MM,SS]= datevec(seconds(timeStep(tim)));
    dateNumber = datenum(Y,M,D,H+initial_time,MM,SS);
        
    plot(lon_epi,lat_epi,'pr','MarkerSize',20,'MarkerFaceColor', 'r');
    lat_ex=[lat_epi-delta_deg lat_epi+delta_deg];
    lon_ex=[lon_epi-delta_deg lon_epi+delta_deg];
    axis([lon_ex lat_ex])
    set(gca, 'CLim', [-0.1,0.1]);
    title(datestr(dateNumber,'HH:MM:SS'));
    colorbar;
    colormap(jet)
    daspect([1 1 1])
    
    
    Frame(tim) = getframe(gcf);
    im = frame2im(Frame(tim)); 
    [imind,cm] = rgb2ind(im,256); 
    
    % Write to the GIF File 
    if tim == 1 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
    clf       
end

% fig =figure;
% movie(fig,Frame,5);


% Create AVI file.
% vidFile = VideoWriter('myPeaks.avi','Motion JPEG AVI');
% vidFile = VideoWriter('newfile.avi','Uncompressed AVI');
% open(vidFile)
% writeVideo(vidFile,Frame)
% close(vidFile)


