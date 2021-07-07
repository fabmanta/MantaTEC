function [tecu,tecuTaped] = filter_tecu_3v(tecuRaw,lowFreq,highFreq,delta)

%apply butterworth filter 1 to 10 mHz to remove daily variation Ionosphere
%(Rolland 2011)  [B,A] = butter(N,Wn)        
Nfilter = 4;                            % 6 filter order
frequency = 1/delta;

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
        
    tecu =nan(size(tecuRaw));    
    for ind_sat = 1:size(tecuRaw,2) % sat loop

        find_gap = sum(isnan(tecuRaw(:,ind_sat)));
        
        if ~(find_gap == length(tecuRaw(:,ind_sat)))
            
            if (find_gap)> 0
                [ ~, block_starts,block_ends,block_lengths] = getBlockIndices(tecuRaw(:,ind_sat),NaN);
                Numb_blocks = length(block_lengths);
            else
                block_starts = 1;
                block_ends = length(tecuRaw(:,ind_sat));
                block_lengths = length(tecuRaw(:,ind_sat));
                Numb_blocks = 1;
            end

            for bn = 1:Numb_blocks
               if block_lengths(bn) <= 24
               tecu(block_starts(bn):block_ends(bn),ind_sat) = nan(block_lengths(bn),1);
               else
                tecu(block_starts(bn):block_ends(bn),ind_sat) = filtfilt(B,A,tecuRaw(block_starts(bn):block_ends(bn),ind_sat));  % get the filtered data   
               end
                    % WITHOUT Signal Processing Toolbox       
                %         [Z,P,G] = myButter(Nfilter,Wnprop,'bandpass');         % find filter coefficients (Alternative)   
                %         [B,A] = zp2tf(Z',P',G);                                %Convert zero-pole-gain filter parameters to transfer function form
                %         tecuTaped = myFiltfilt(B,A,tecuRawTaped,[]);
            end
        else
           tecu(:,ind_sat) = tecuRaw(:,ind_sat);  % get the filtered data   

        end

    end
end

