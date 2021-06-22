function [ blocks, block_starts, block_ends, block_lengths ] = getBlockIndices( input, gap_value )
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

blocks = find(~isnan(input));

if any(blocks)
    %% Define the window over which there are gaps
    clear block_starts block_ends; % make sure gap_starts and gap_ends are clear
    block_starts = blocks(1); % initialize gap_starts as the index for the first gap
    block_ends = []; % initialize gap_ends as blank
    for n = 1:length(blocks)-1
        if(blocks(n+1)-blocks(n)) > 1
            block_starts = [block_starts; blocks(n+1)];
            block_ends = [block_ends; blocks(n)];
        end
    end
    block_ends = [block_ends; blocks(end)]; % finalize the gap_ends with the index for the last gap
    block_lengths = block_ends - block_starts + 1;
else
    blocks = 0;
    block_starts = NaN; 
    block_ends = NaN;
    block_lengths = 0;
    
end

end
