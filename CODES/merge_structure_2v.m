function data = merge_structure_2v(t,sip,tec,sv)

% MAKE_DATA	Merge sr, sf, sip, and lgu structures into a single structure 'data'
%                 sf = filtered IEC
%                 sr = unfiltered IEC
%                 sl = LG corrected for ambiguity bias but not for interfequency bias
%                 sip = 
%		Structure 'data' contains, for each prn:
%		   time_iec iec_filt iec_raw time_sip lat_sip lon_sip elev azim dist_epi azim_epi mag_los emf lgu
%		sv is a vector with the list of prn's to consider
%
%		data = make_data(sr,sf,sl,sip,sv);
%

% message
disp(['Merging iec and sip into single data structure']);

% initialize structure
data = [];

% for each prn
for i=1:length(sv)

  % extract data
  eval(['tmp_tec = tec.' (sv(i,:)) ';']);
  eval(['tmp_sip = sip.' (sv(i,:)) ';']);

  % extract time vector and round to the second
  te = round(t*3600);
  ts = round(tmp_sip(:,1)*3600);

  % find sip entries matching filtered iec (end of iec may not have sip values?)
  I = []; 
  for j=1:length(te)
    I = [I find(ts==te(j))];
  end

  % rewrite tmp_fil
  tmp_tec = tmp_tec(1:length(I),:);

  % extract arrray from matching sip
  match_sip = tmp_sip(I,:);


  % fill up data structure
  tec_sip = [match_sip tmp_tec];
  
  data = setfield(data,sv(i,:),tec_sip);
  
end
