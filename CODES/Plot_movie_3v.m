%% PLOT IPP video
% Author :Fabio Manta 
% Last update: March-2021

id_staz = (1:length(staz));
SatNum = 1:size(name_sat_srt,1); % all sats;
% SatNum = [2 5 6 11 16 18 23 24 27 39 40 40 50];
filename = 'movie_3.gif';
initial_time = 19.25; % in hours
dt= 30;
timeStep = (0:dt:3600*0.75);
delta_deg =18;
lat_ex=[lat_epi-delta_deg lat_epi+delta_deg];
lon_ex=[lon_epi-delta_deg lon_epi+delta_deg];
type = 2; % 1 = diffTEC 2=vtec
radius = 6371;
PPI = 300;

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
            
            if (isfield(eval(staz(id_staz(jj),1:end-4)),sprintf('%s',name_sat_srt(SatNum(ii),:))  ))
         % plot all station tecu
            tmpTIME = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,1)',name_sat_srt(SatNum(ii),:)))); 
            tmpLAT = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,2)',name_sat_srt(SatNum(ii),:))));  
            tmpLON = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,3)',name_sat_srt(SatNum(ii),:))));
            tmpELE = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,4)',name_sat_srt(SatNum(ii),:))));

            if (type == 1)
                tmpTecu = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,7)',name_sat_srt(SatNum(ii),:))));                       
                tmpTecu = diff(tmpTecu);
                
                tmpLON = tmpLON(1:length(tmpTecu));
                tmpLAT = tmpLAT(1:length(tmpTecu));
                tmpELE = tmpELE(1:length(tmpTecu));

                alpha = (radius.*cosd(tmpELE))./(radius+PPI);
                obliquity = sqrt(1-(alpha.^2));
                tmpTecu = tmpTecu.*obliquity;
                A = 0.1;
            elseif (type == 2)
                tmpTecu = eval(strcat(staz(id_staz(jj),1:end-4),sprintf('.%s(:,8)',name_sat_srt(SatNum(ii),:))));                       
                tmpTecu = signal_attenuation(tmpTecu);
                A = 0.07;
            end
               
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
            end% if the satellite exist 
        end
    end
    
    dateNumber = datenum(0,0,0,reference_time,0,0);
    [xEvent, yEvent]= m_ll2xy(lon_epi,lat_epi);
    plot(xEvent,yEvent,'pr','MarkerSize',20,'MarkerFaceColor', 'r','MarkerEdgeColor','k');

    set(gca, 'CLim', [-A,A]);
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
