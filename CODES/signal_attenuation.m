function tec_smooth = signal_attenuation(tmpTEC)
% figure;plot(tmpTEC);
n_points = 10;
% copy the last XX point of the time seires and reverse them
% then normalize in between O and 1
% multiply it fo the corresponding point index in the timeseries


    % find gaps index
    [gaps, gap_starts, gap_ends, gap_lengths] = getGapIndices(tmpTEC,NaN);
    att_vect = linspace(0,1,n_points+1)';
    att_vect = linspace(0,n_points,n_points+1)';
    Y = exp(att_vect/2);
    att_vect = (Y - min(Y)) / ( max(Y) - min(Y) );
    
if ~isnan(gap_starts)
        
        for gg = 1:size(gap_lengths,1)
            
            if gap_starts(gg) == 1
                ini =gap_starts(gg)+gap_lengths(gg);         
                if size(tmpTEC(ini:end ),1)> n_points+1 % if the data points after the gap are more then the length of the gap then...
                    tmpTEC(ini:ini+n_points) = tmpTEC(ini:ini+n_points).*att_vect;
                end
            else
                ini =gap_starts(gg)-n_points;         
                if size(tmpTEC(ini:end ),1)> n_points+1 % if the data points after the gap are more then the length of the gap then...
                    
                    tmpTEC(gap_starts(gg)-1-n_points:gap_starts(gg)-1) = tmpTEC(gap_starts(gg)-1-n_points:gap_starts(gg)-1).*flip(att_vect);
                    if gap_ends(gg) < length(tmpTEC)
                        tmpTEC(gap_ends(gg)+1:gap_ends(gg)+1+n_points) = tmpTEC(gap_ends(gg)+1:gap_ends(gg)+1+n_points).*(att_vect);
                    end
                    
                end
            end
            
        end
    
end
tec_smooth = tmpTEC;

