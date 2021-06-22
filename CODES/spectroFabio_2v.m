function spectroFabio_2v(TIME,filtTECU,GNSSname,n_satel,EventTimeHours,fs,wind,nfft,dist,oneHourBefore,twoHourAfter,Tmin)
   
        GNSS = upper(GNSSname);
        filtTECU_ungapped = fillmissing(filtTECU,'previous');
        filtTECU_ungapped(isnan(filtTECU_ungapped)) = 0;

        Hours = floor(EventTimeHours);
        fract = (EventTimeHours-Hours);
        Minutes = floor(fract*60);
        Seconds = ((fract*60)- Minutes)*60;
        EventTimeSecs = Hours*3600+ Minutes*60 + Seconds;
        TminSec = Tmin*3600;

       
    figure('units','normalized','outerposition',[0 0 1 1])
    set(gcf,'color','w');
    
    subplot(211);
        set(gca,'FontSize',16)
        hold all
        plot(TIME,filtTECU,'k','LineWidth',3) 
        ylimit = 0.5;
        box on
        plot([EventTimeHours EventTimeHours],[-ylimit,ylimit],'r','LineWidth',1.5)

        plot([EventTimeHours+(8/60) EventTimeHours+(8/60)],[-ylimit, ylimit],'--','LineWidth',1, 'Color',[0.1 0.1 0.1])
        axis([oneHourBefore twoHourAfter -ylimit ylimit])
        set(gca,'XTick',(oneHourBefore:0.5:twoHourAfter))
        annotation('textbox',[0.68,0.83,0.1,0.1],'String',(sprintf('%s-PRN%d @ %4.0fkm',GNSS(1:4),n_satel,dist)),'FontSize',18,'EdgeColor','none','FontWeight','Bold')
        xlabel('Time (UTC)','FontSize', 18,'FontWeight','Bold');
        ylabel('Amplitude (TECU)','FontSize', 18,'FontWeight','Bold');
        
        
  sub = subplot(212);
        set(gca,'FontSize',16)
        hold all
        set(sub,'Position',[0.1300 0.1450 0.7750 0.3412]);
        [S,F,T,P] = spectrogram(filtTECU_ungapped,wind,wind-1,nfft,fs,'yaxis'); 
        P2 = abs(P+eps)./400;
        pdb = 10*log10(P2/1e-12); %convert power measurements in decibels

        pcolor(T+TminSec,F*1.e+3,pdb)
        shading flat
        set(gca, 'CLim', [80 110]); %default 90 110
%         set(gca, 'CLim', [70 87]);
        set(gca, 'YTick',[0 2 5 10]);
        ylim([0 10])
        xlim([oneHourBefore*3600, twoHourAfter*3600])

        set(gca,'xtick',[])
        colormap(hot)
        c =colorbar('southoutside');
        set(c,'position',[0.1304 0.095 0.7748 0.0347],'FontSize', 18)
        
        plot([EventTimeSecs EventTimeSecs],[0, 10],'w','LineWidth',1.5)
        plot([EventTimeSecs+480 EventTimeSecs+480],[0, 10],'--w','LineWidth',1) %8min after
        
        v =axis;
        plot([v(1) v(2)],[2, 2],'--w','LineWidth',1)
        ylabel('Frequency (mHz)','FontSize', 18,'FontWeight','Bold');
end

