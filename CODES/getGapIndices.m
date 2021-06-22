function [ gaps, gap_starts, gap_ends, gap_lengths ] = getGapIndices( input, gap_value )
%{
Parses through a vector of data and returns the indices for gaps in the
data. Gaps are defined by a specific value as well as any NaN values.
Output fields are:
    gaps  :  a vector of the indices that are gaps
    gap_starts  :  a vector with the index for the beginning of each gap
                   window
    gap_ends     :  a vector with the index for the end of each gap window
    gap_lengths  :  a vector with the length of each gap window
                    (gap_ends - gap_starts + 1)

%}

%% Create a vector of gap indices
gaps = find(input==gap_value); % creates a vector with gaps (gaps defined as zero)
gaps = [gaps; find(isnan(input))]; % adds NaN values to the gap vector
gaps = sort(gaps); % puts the zero and NaN gap points in order

if any(gaps)
    %% Define the window over which there are gaps
    clear gap_starts gap_ends; % make sure gap_starts and gap_ends are clear
    gap_starts = gaps(1); % initialize gap_starts as the index for the first gap
    gap_ends = []; % initialize gap_ends as blank
    for n = 1:length(gaps)-1
        if(gaps(n+1)-gaps(n)) > 1
            gap_starts = [gap_starts; gaps(n+1)];
            gap_ends = [gap_ends; gaps(n)];
        end
    end
    gap_ends = [gap_ends; gaps(end)]; % finalize the gap_ends with the index for the last gap
    gap_lengths = gap_ends - gap_starts + 1;
else
    gaps = 0;
    gap_starts = NaN; 
    gap_ends = NaN;
    gap_lengths = 0;
    
end

end
