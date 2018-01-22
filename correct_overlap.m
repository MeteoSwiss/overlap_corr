function chm=correct_overlap(chm,station,serlom,folder_overlap_model)
% chm=correct_overlap(chm,serlom,serlom,folder_overlap_model)
% Corrects the variable beta_raw for the overlap artefacts.
% The overlap function is modeled according to the internal temperature
% More info: http://www.atmos-meas-tech-discuss.net/amt-2016-30/
% hem/poy Meteoswiss 2016

if nargin < 4
    if isunix()
        folder_overlap_model = '/data/pay/PBL4EMPA/overlap_correction/';
    else
        folder_overlap_model = 'M:\pay-data\data\pay\PBL4EMPA\overlap_correction\';
    end
end

%% Load file
file_overlap = [folder_overlap_model '/' 'Overlap_correction_model_' station '_' serlom  '.nc'];
if exist (file_overlap,'file')==0
    warning(['No overlap correction file found: ' file_overlap])
    return
else
    disp([ 'Read : ' file_overlap])
end
a=ncread(file_overlap,'a');
b=ncread(file_overlap,'b');
overlap_ref=ncread(file_overlap,'overlap_ref');

%% Correct overlap
disp('Correct overlap')
RCS_raw=chm.beta_raw;
ov_rec_all = NaN(size(RCS_raw));
%     for i=1:find(chm.range<=1200,1,'last')
for i=1:length(chm.range)
    rel_diff = polyval([a(i),b(i)],chm.temp_int-273.15);
    ov_rec_all(i,:) = 1./ (rel_diff/100/overlap_ref(i) + 1./overlap_ref(i));
end
%     ov_rec_all(find(chm.range<=1200,1,'last')+1:end,:) = repmat(overlap_ref(find(chm.range<=1200,1,'last')+1:end),1,size(RCS_raw,2));

chm.beta_raw_ovcorr = chm.beta_raw./ov_rec_all.*repmat(overlap_ref,1,size(RCS_raw,2));
chm.beta_raw_ovcorr_0_05ov = chm.beta_raw_ovcorr;
chm.beta_raw_ovcorr_0_05ov(overlap_ref<0.05,:) = NaN;
chm.ovcorr = ov_rec_all;

