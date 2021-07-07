%% Cycle-slip detection and fix
% Author :Fabio Manta 
% Last update: July-2021
function [x2] = cycle_slip(t,x1,maxnumchg)

%% Example
% clear all
% t = 0 : 0.1 : 10*pi;
% x1 = sin(t);
% x1(1,131:315) = x1(1,131:315) + 5;
% x1(1,170:315) = x1(1,170:315) - 10;
% x1 = x1+ randn(size(t));
% 
% 
% % Find the gaps
% maxnumchg = 4; % max number of change
% [ipoint,residual] = findchangepts(x1,...
%    'MaxNumChanges',maxnumchg);
% 
% if any(ipoint)
%     % Define the window over which there are gaps
%     clear block_starts block_ends; % make sure gap_starts and gap_ends are clear
%     block_starts = 1; % initialize gap_starts at the beginning of the time series
%     block_ends = []; % initialize gap_ends as blank
%     for n = 1:length(ipoint)
% %         if(ipoint(n+1)-ipoint(n)) > 1
%             block_starts = [block_starts; ipoint(n)];
%             block_ends = [block_ends; ipoint(n)-1];
% %         end
%     end
%     block_ends = [block_ends; length(x1)]; % finalize the gap_ends with the index for the last gap
%     block_lengths = block_ends - block_starts + 1;
% else
%     ipoint = 0;
%     block_starts = NaN; 
%     block_ends = NaN;
%     block_lengths = 0;
%     
% end
% 
% for ii = 1 : length(block_lengths)
%     m(ii) = mean(x1(block_starts(ii):block_ends(ii)));
% end
% 
% x2 = zeros(size(x1));
% for ii = 1 : length(block_lengths)
%     x2(block_starts(ii):block_ends(ii)) = x1(block_starts(ii):block_ends(ii)) - m(ii);
% end
% 
% figure
% subplot(311)
% plot(t,x1,'b'); title('Original Data with Jumps')
% 
% subplot(312)
% plot(t,x1,'b')
% hold on
% for ii = 1 : length(block_lengths)
%     plot([t(block_starts(ii)) t(block_ends(ii))],[m(ii) m(ii)],'k')
% end
% title('Original Data with Change Points')
% 
% subplot(313)
% plot(t,x1,t,x2)
% title('Data without Jumps')


%%

% Find the gaps
[ipoint,residual] = findchangepts(x1,...
   'MaxNumChanges',maxnumchg);

if any(ipoint)
    % Define the window over which there are gaps
    clear block_starts block_ends; % make sure gap_starts and gap_ends are clear
    block_starts = 1; % initialize gap_starts at the beginning of the time series
    block_ends = []; % initialize gap_ends as blank
    for n = 1:length(ipoint)
%         if(ipoint(n+1)-ipoint(n)) > 1
            block_starts = [block_starts; ipoint(n)];
            block_ends = [block_ends; ipoint(n)-1];
%         end
    end
    block_ends = [block_ends; length(x1)]; % finalize the gap_ends with the index for the last gap
    block_lengths = block_ends - block_starts + 1;
else
    ipoint = 0;
    block_starts = NaN; 
    block_ends = NaN;
    block_lengths = 0;
    
end

for ii = 1 : length(block_lengths)
    m(ii) = mean(x1(block_starts(ii):block_ends(ii)));
end

x2 = zeros(size(x1));
for ii = 1 : length(block_lengths)
    x2(block_starts(ii):block_ends(ii)) = x1(block_starts(ii):block_ends(ii)) - m(ii);
end

