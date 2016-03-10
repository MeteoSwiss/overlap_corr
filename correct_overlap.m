function chm=correct_overlap(chm,serlom,folder_overlap_model)
% chm=correct_overlap(chm,serlom,folder_overlap_model)
% Corrects the variable beta_raw for the overlap artefacts.
% The overlap function is modeled according to the internal temperature

file_overlap = [folder_overlap_model '/' 'Overlap_correction_model_' serlom  '.nc'];

a=ncread(file_overlap,'a');
b=ncread(file_overlap,'b');
overlap_ref=ncread(file_overlap,'overlap_ref');

ov_rec_all = NaN(size(RCS_raw));
%     for i=1:find(chm.range<=1200,1,'last')
for i=1:length(chm.range)
    rel_diff = polyval([a(i),b(i)],chm.temp_int-273.15);
    ov_rec_all(i,:) = 1./ (rel_diff/100/overlap_ref(i) + 1./overlap_ref(i));
end
%     ov_rec_all(find(chm.range<=1200,1,'last')+1:end,:) = repmat(overlap_ref(find(chm.range<=1200,1,'last')+1:end),1,size(RCS_raw,2));

chm.beta_raw_ovcorr = chm.beta_raw/ov_rec_all.*repmat(overlap_ref,1,size(RCS_raw,2));
chm.ovcorr = ov_rec_all;