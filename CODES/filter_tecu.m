function [tecu,tecuTaped] = filter_tecu(tecuRaw,lowFreq,highFreq,delta)

        %apply butterworth filter 1 to 10 mHz to remove daily variation Ionosphere
        %(Rolland 2011)  [B,A] = butter(N,Wn)
        
        for ind_sat = 1:size(tecuRaw,2) % sat loop
            tmp_tecuRaw = [];
            tmp_tecuRaw = tecuRaw(:,ind_sat);
            find_gap = sum(isnan(tmp_tecuRaw));
            if (find_gap)> 0
                [~, gap_starts, ~, gap_lengths] = getGapIndices(tmp_tecuRaw,NaN);
                tecuRawTaped(:,ind_sat) = fillgapstaper(tmp_tecuRaw);
            else
                tecuRawTaped(:,ind_sat) = tmp_tecuRaw;
            end
        end
        
        frequency = 1/delta;
        Nfilter = 4;                            % 6 filter order
 
        % WITH Signal Processing Toolbox   
        if (delta == 120)
            Wn = 0.004;                            %best used for 2004 EQ SuGAr data
            Wnprop = Wn./(frequency/2);               % Wn must be smaller than 1, do the proportion (1 = frequency/2)
            [B,A] = butter(Nfilter,Wnprop,'high');  % find filter coefficients (use "high" only for Sumatra 2004)
        else
            Wn = [lowFreq,highFreq];
            Wnprop = Wn./(frequency/2);               % Wn must be smaller than 1, do the proportion (1 = frequency/2)
            [B,A] = butter(Nfilter,Wnprop);           % find filter coefficients (use "high" only for Sumatra 2004)
        end
        tecuTaped = filtfilt(B,A,tecuRawTaped);  % get the filtered data   

        % WITHOUT Signal Processing Toolbox       
%         [Z,P,G] = myButter(Nfilter,Wnprop,'bandpass');         % find filter coefficients (Alternative)   
%         [B,A] = zp2tf(Z',P',G);                                %Convert zero-pole-gain filter parameters to transfer function form
%         tecuTaped = myFiltfilt(B,A,tecuRawTaped,[]);

        tecu = tecuTaped;
        
        
        for ind_sat = 1:size(tecuRaw,2)
            tmp_tecuRaw = [];
            tmp_tecuRaw = tecuRaw(:,ind_sat);
            find_gap = sum(isnan(tmp_tecuRaw));
            
            if (find_gap)> 0
                [~, gap_starts, ~, gap_lengths] = getGapIndices(tmp_tecuRaw,NaN);
                for gapp =1:length(gap_starts)
                    tecu(gap_starts(gapp):gap_starts(gapp)+gap_lengths(gapp)-1,ind_sat) = NaN; 
                end
            end
        end

end

