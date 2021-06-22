%% PLOT IPP video
% Author :Fabio Manta 
% Last update: April-2021


function Plot_movie_4v(proj_type,lat_ex,lon_ex,lat_epi,lon_epi,id_staz,name_sat_srt,SatNum,initial_time,dt,timeStep,radius,PPI,type,filename,outputDirectory,plotsDirectory)
    
global staz; % list of GNSS stations
    figure;

    for tim = 1:size(timeStep,2)

        m_proj(proj_type,'longtitudes',lon_ex,'latitudes',lat_ex);
        hold on
        m_coast('color',[0 0 0]); 
        shading flat;
        m_grid('linest','none','box','fancy','tickdir','in','fontsize',16);
        hold on

        for jj = 1:size(id_staz,2)
            
            fileExist =(staz(id_staz(jj),1:8));
            load(fullfile(outputDirectory,fileExist));

            for ii = 1:length(SatNum)  
                hold on    

                if isfield(eval(fileExist),sprintf('%s',name_sat_srt(SatNum(ii),:)))
                    
             % plot all station tecu
                tmpTIME = eval(strcat(staz(id_staz(jj),1:8),sprintf('.%s(:,1)',name_sat_srt(SatNum(ii),:)))); 
                tmpLAT = eval(strcat(staz(id_staz(jj),1:8),sprintf('.%s(:,2)',name_sat_srt(SatNum(ii),:))));  
                tmpLON = eval(strcat(staz(id_staz(jj),1:8),sprintf('.%s(:,3)',name_sat_srt(SatNum(ii),:))));
                tmpELE = eval(strcat(staz(id_staz(jj),1:8),sprintf('.%s(:,4)',name_sat_srt(SatNum(ii),:))));

                    if (type == 1)
                        tmpTecu = eval(strcat(staz(id_staz(jj),1:8),sprintf('.%s(:,7)',name_sat_srt(SatNum(ii),:))));                       
                        tmpTecu = diff(tmpTecu);

                        tmpLON = tmpLON(1:length(tmpTecu));
                        tmpLAT = tmpLAT(1:length(tmpTecu));
                        tmpELE = tmpELE(1:length(tmpTecu));

                        alpha = (radius.*cosd(tmpELE))./(radius+PPI);
                        obliquity = sqrt(1-(alpha.^2));
                        tmpTecu = tmpTecu.*obliquity;
                        A = 0.1;
                    elseif (type == 2)
                        tmpTecu = eval(strcat(staz(id_staz(jj),1:8),sprintf('.%s(:,8)',name_sat_srt(SatNum(ii),:))));                       
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

end