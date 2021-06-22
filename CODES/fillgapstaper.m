function [ output ] = fillgapstaper( data, varargin )
%{
% Fills data gaps by applying a taper. The taper is the product of non-gap
% values before a given gap and after a given gap. The first step of this
% script is to smooth over any data gaps that are only one sample long. The
% default gap value is NaN. An additional gap value can be passed to the
% function by including it as an input variable after the vector.
% 
% Examples:
% >> filledvector1 = fillgapstaper( vector1 ); % fills all NaN values
% 
% >> filledvector2 = fillgapstaper( vector2, 0); fills all NaN and 0 values
%}

%% Parse inputs for additional gap values, if there are any

% If the number of input arguments is greater than 1 (the number of 
% explicitly declared variables), define gap_value as the first variable
% input. Otherwise, set the gap_value to NaN as a default.

if nargin > 1
    gap_value = varargin{1};
else
    gap_value = NaN;
end


%% Smooth over 1 sample gaps

% First get gap indices
[~, gap_starts, ~, gap_lengths] = getGapIndices(data,gap_value);

% Now replace 1 sample gaps with averaged point
for n = 1:length(gap_starts)
    
    if(gap_lengths(n)==1)
        index = gap_starts(n);
        pre_index = index-1;
        post_index = index+1;
        
        if (pre_index == 0)
            data(index) =  data(post_index);   
        elseif (post_index >  length(data))
            data(index) = data(pre_index);   
        else
            data(index) = (data(pre_index)+data(post_index))/2;
        end
    end
    
end


%% Fill Gaps with a taper
  find_gap = sum(isnan(data));
        if (find_gap)> 0

            % Get the location and indices of the gaps again
            [~, gap_starts, gap_ends, ~] = getGapIndices(data,gap_value);

            % Create and apply the tapers
            for n = 1:length(gap_starts) % for each gap

                front = gap_starts(n); % get the beginning of the gap
                back = gap_ends(n); % get the end of the gap
                gap_length = back-front+1; % compute the gap length

                fsegment = collapseData(data(1:front-1),0); % collapses all data before the gap
                
                if isempty(fsegment)
                    fsegment = zeros(1,gap_length);
                else
                    if length(fsegment) < gap_length

                        Numb_times = floor(gap_length/length(fsegment));
                        patching = fliplr(fsegment');
                        fsegment = zeros(1,gap_length);        
                        patching = repmat(patching,1,Numb_times);
                        fsegment(1:length(patching)) = patching;

                    else
                        fsegment = fsegment(end-gap_length+1:end); % takes the last w samples before the gap, where w is the gap_length
                        fsegment = fliplr(fsegment'); % reverse the order of the segment
                    end

                    fsegment = fsegment.*linspace(1,0,length(fsegment')); % multiplies by a line from 1 to 0
                end
            %% 
                patching = [];
                bsegment = collapseData(data(back+1:end),0); % collapses all data before the gap
                
                if isempty(bsegment)
                    bsegment = zeros(1,gap_length);  
                else
                    if  length(bsegment) < gap_length
                        Numb_times = floor(gap_length/length(bsegment));
                        patching = fliplr(bsegment');
                        bsegment = zeros(1,gap_length);        
                        patching = repmat(patching,1,Numb_times);
                        bsegment(end-length(patching)+1:end) = patching;
                    else
                        bsegment = bsegment(1:gap_length); % takes the last w samples before the gap, where w is the gap_length
                        bsegment = fliplr(bsegment'); % reverse the order of the segment
                    end
                    bsegment = bsegment.*linspace(0,1,length(bsegment')); % multiplies by a line from 0 to 1
                end
                
                
                segment = fsegment + bsegment;
                data(front:back) = segment;

            end
            
        end
%% Output the data

output = data;

end
%% :)