%% PLOT IPP video
% Author :Fabio Manta 
% Last update: Dec-2020

id_staz = (1:length(staz));
SatNum = [18 25 26 29 31];
filename = 'movie_2.gif';
initial_time = 17; % in hours
dt= 30;
timeStep = (0:dt:3600*1);
delta_deg =6;
lat_ex=[lat_epi-delta_deg lat_epi+delta_deg];
lon_ex=[lon_epi-delta_deg lon_epi+delta_deg];


figure;

for tim = 1:size(timeStep,2)

%     m_proj('lambert','longtitudes',Lon_limit,'latitudes',Lat_limit);
    m_proj('Mercator','longtitudes',lon_ex,'latitudes',lat_ex);
    hold on
    m_coast('color',[0 0 0]); 
    shading flat;
    m_grid('linest','none','box','fancy','tickdir','in','fontsize',16);


        
    hold on
    
    for jj = 1:size(id_staz,2)
         
        for ii = 1:length(SatNum)  
            hold on    
         % plot all station tecu
            tmpTIME = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,1)',SatNum(ii)))); 
            tmpLAT = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,2)',SatNum(ii))));  
            tmpLON = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,3)',SatNum(ii))));
            
            tmpTecu = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.PRN%d(:,8)',SatNum(ii))));                       
    
               
            sec2hours = timeStep(tim)/3600;
            reference_time = initial_time+sec2hours;
            startingPosition = knnsearch(tmpTIME, reference_time);
            istime = abs(tmpTIME(startingPosition)-reference_time)< dt;

                    if (istime)== 1

                        if startingPosition+tim-1<=(length(tmpLON))
                                       
                            [xipp, yipp]= m_ll2xy(tmpLON(startingPosition),tmpLAT(startingPosition));
                            cplot(xipp, yipp...
                            ,tmpTecu(startingPosition),'o','markerfacecolor','flat','MarkerEdgeColor','none','MarkerSize',8);
                        
                        end
                     
                        
                    end

        end
    end
    
    dateNumber = datenum(0,0,0,reference_time,0,0);
    [xEvent, yEvent]= m_ll2xy(lon_epi,lat_epi);
    plot(xEvent,yEvent,'pr','MarkerSize',20,'MarkerFaceColor', 'r','MarkerEdgeColor','k');

    set(gca, 'CLim', [-0.09,0.09]);
    title(datestr(dateNumber,'HH:MM:SS'));
    colorbar;
    colormap(jet)
       
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


close all
