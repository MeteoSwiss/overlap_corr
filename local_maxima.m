function [val,ind] = local_maxima(v,idx)

idx = 1:length(v);
idx = idx(idx>2 & idx<length(v)-1);

val = NaN(0,2);
ind = NaN(0,2);

for i=1:length(idx) 
    if v(idx(i))>v(idx(i)-1) && v(idx(i))>v(idx(i)+1) ...
    && v(idx(i))>v(idx(i)+2) && v(idx(i))>v(idx(i)-2)
            val(end+1) = v(idx(i)); %value
            ind(end+1) = idx(i); %index
    end
end
% for i=1:length(idx) 
%     if v(idx(i))>=nanmean([v(idx(i)-1) v(idx(i)-2)]) && ...
%        v(idx(i))>=nanmean([v(idx(i)+1) v(idx(i)+2)])
%             val(end+1) = v(idx(i)); %value
%             ind(end+1) = idx(i); %index
%     end
% end
end