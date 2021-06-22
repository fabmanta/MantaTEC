function vtec = vertical_tecu_calculator_2v(iec,SIP_mat,PPI,sv)

radius = 6371;
PPI =PPI/1000;
vtec = nan(size(iec));

for j= 1:size(iec,2)
    
    tf = isstruct(SIP_mat);
    if tf == 1
        eval(['elev = SIP_mat.PRN' num2str(sv(j)) '(:,4);']);
    else
        elev = SIP_mat;
    end
    
    a = (radius.*cosd(elev))./(radius+PPI);
    obliquity = sqrt(1-(a.^2));
    vtec(:,j) = iec(:,j).*obliquity;
end

end