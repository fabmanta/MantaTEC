function [tecu] = filter_tecu_2v(tecuRaw,lowFreq,highFreq,delta)

        %apply butterworth filter 1 to 10 mHz to remove daily variation Ionosphere
        %(Rolland 2011)  [B,A] = butter(N,Wn)
        
        
        frequency = 1/delta;
        Nfilter = 4;                            % 6 filter order
 
        % WITH Signal Processing Toolbox   
        if (delta == 120)
            Wn = 0.001;                            %best used for 2004 EQ SuGAr data
            Wnprop = Wn./(frequency/2);               % Wn must be smaller than 1, do the proportion (1 = frequency/2)
            [B,A] = butter(Nfilter,Wnprop,'high');  % find filter coefficients (use "high" only for Sumatra 2004)
        else
            Wn = [lowFreq,highFreq];
            Wnprop = Wn./(frequency/2);               % Wn must be smaller than 1, do the proportion (1 = frequency/2)
            [B,A] = butter(Nfilter,Wnprop);           % find filter coefficients (use "high" only for Sumatra 2004)
        end
        tecuRaw(~isnan(tecuRaw)) = filtfilt(B,A,tecuRaw(~isnan(tecuRaw)));  % get the filtered data   

        % WITHOUT Signal Processing Toolbox       
%         [Z,P,G] = myButter(Nfilter,Wnprop,'bandpass');         % find filter coefficients (Alternative)   
%         [B,A] = zp2tf(Z',P',G);                                %Convert zero-pole-gain filter parameters to transfer function form
%         tecuTaped = myFiltfilt(B,A,tecuRawTaped,[]);

        tecu = tecuRaw;
        
      

end

