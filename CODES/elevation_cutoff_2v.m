function [SIP_mat,TEC_mat] = elevation_cutoff_2v(SIP_mat,TEC_mat,ele_cutoff,name_sv,sv)


% t(index_cutoff) = [];
for sat_n = 1 : size(sv,1)
    index_cutoff =[];
    if strcmp(name_sv(sat_n,1),'G')
        control = eval(sprintf('SIP_mat.G%.2d(:,4)',sv(sat_n)));
        index_cutoff = find( control < ele_cutoff); % find elements with elevation lower than the cutoff 
        eval([sprintf('SIP_mat.G%.2d',sv(sat_n)) '(index_cutoff,2:end) =NaN;']);
        eval([sprintf('TEC_mat.G%.2d',sv(sat_n)) '(index_cutoff,2:end) =NaN;']);
    elseif strcmp(name_sv(sat_n,1),'R')
        control = eval(sprintf('SIP_mat.R%.2d(:,4)',sv(sat_n)));
        index_cutoff = find( control < ele_cutoff); % find elements with elevation lower than the cutoff 
        eval([sprintf('SIP_mat.R%.2d',sv(sat_n)) '(index_cutoff,2:end) =NaN;']);
        eval([sprintf('TEC_mat.R%.2d',sv(sat_n)) '(index_cutoff,2:end) =NaN;']);        
    end
end