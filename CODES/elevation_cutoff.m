function [SIP_mat,TEC_mat] = elevation_cutoff(SIP_mat,TEC_mat,ele_cutoff,sv)


% t(index_cutoff) = [];
for sat_n = 1 : size(sv,1)
    index_cutoff =[];
    control = eval(sprintf('SIP_mat.PRN%d(:,4)',sv(sat_n)));
    index_cutoff = find( control < ele_cutoff); % find elements with elevation lower than the cutoff 

    eval(['SIP_mat.PRN' num2str(sv(sat_n)) '(index_cutoff,2:end) =NaN;']);
    eval(['TEC_mat.PRN' num2str(sv(sat_n)) '(index_cutoff,:) =NaN;']);
end