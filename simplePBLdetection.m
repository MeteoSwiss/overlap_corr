function [gradients,xtime] = simplePBLdetection(chm,overlap_ref,overlap_to_use,do_plot)
% Calculate PBL with simple code
% gradients(1,:) = lowest gradient
% gradients(4:-1;2,:) = 3 strongest gradients
% gradients(5,:) = lowest cloud base

if nargin==3
   do_plot = 0; 
end

chm.signal_raw_per_pulse = chm.beta_raw .* repmat(overlap_ref,1,size(chm.beta_raw,2)) ./ repmat(chm.range.^2,1,size(chm.beta_raw,2)) .* (chm.scaling*repmat(chm.p_calc',size(chm.beta_raw,1),1)) + repmat(chm.base',size(chm.beta_raw,1),1);

p = polyfit(sqrt(chm.base),chm.stddev,1);
chm.noise_level = p(1) * sqrt(chm.signal_raw_per_pulse);
chm.SNR = abs((chm.signal_raw_per_pulse - repmat(chm.base',size(chm.beta_raw,1),1))./chm.noise_level);

A = sparse(diag(overlap_ref./chm.range.^2));
alpha = 2.5e-16;
chm.beta_raw_2 = (A'*A+alpha*speye(size(A)))\(A'*(chm.beta_raw .* repmat(overlap_ref,1,size(chm.beta_raw,2)) ./ repmat(chm.range.^2,1,size(chm.beta_raw,2))));
chm.noise_level_2 = (A'*A+alpha*speye(size(A)))\(A'*(chm.noise_level ./ (chm.scaling*repmat(chm.p_calc',size(chm.beta_raw,1),1))));
chm.SNR_2 = abs(chm.beta_raw_2 ./ chm.noise_level_2);

if numel(overlap_to_use) == length(chm.time)*length(chm.range)
    RCS = chm.beta_raw_2./overlap_to_use.*repmat(overlap_ref,1,size(chm.beta_raw_2,2));
else
    RCS = chm.beta_raw_2./repmat(overlap_to_use,1,size(chm.beta_raw_2,2)).*repmat(overlap_ref,1,size(chm.beta_raw_2,2));
end
SNR = chm.SNR_2;

SNRlim = 1.96;
maskSNR = (SNR>SNRlim).*(RCS>0);
for j=1:3
    ext = [ones(1,length(chm.time)+2);ones(length(chm.range),1),maskSNR,ones(length(chm.range),1);ones(1,length(chm.time)+2)];
    cv = conv2(ext,ones(3,3),'valid');
    maskSNR(cv~=9) = 0;
end
for j=1:15
    cv = conv2(maskSNR,ones(3,3),'same');
    maskSNR(cv>0) = 1;
end
maskUC = maskSNR;
for j=1:length(chm.time)
   null_tmp = cumsum(~maskSNR(:,j));
   one_tmp = cumsum(maskSNR(:,j));
   tmp = one_tmp - null_tmp;
   dtmp = conv(tmp,[1;0;-1],'same');
   d2tmp = conv(dtmp,[1;0;-1],'same');
   ind_p = find((dtmp==0).*(d2tmp<0).*(chm.range>600),1,'first');
   if ~isempty(ind_p)
      maskUC(ind_p+2:end,j) = 0; 
   end
end

xtime = 0:10/60/24:1;
maxR = 2500;
minR = chm.range(find(overlap_ref>=0.1,1,'first'));
is_bad_weather = chm.sci~=0;
gradients = NaN(5,length(xtime));

if do_plot
    figure('NextPlot','ReplaceChildren');
end

for j=1:length(xtime)-1
   indt = chm.time-floor(chm.time(1))>=xtime(j) & chm.time-floor(chm.time(1))<xtime(j+1);
   
   if any(is_bad_weather(indt))
       continue;
   end
   
   disp(datestr(floor(chm.time(1))+xtime(j+1)));
   
   profile = nanmean(log10(abs(RCS(:,indt))),2);
   
   cbh = chm.cbh(1,indt)-chm.cho;cbh(cbh<0) = NaN;
   maxReffective = min(maxR,nanmin(cbh));
   
   gradients(5,j+1) = nanmin(cbh);
   
   rSNR = chm.range(find(chm.range>=600 & nanmin(maskUC(:,indt),[],2)<1,1,'first'));
   maxReffective = min(maxReffective,rSNR);
   
   if maxReffective<minR
       continue;
   end
   
   
   [~,grad_profile,~] = sgolay_smooth(chm.range(1:find(chm.range<=maxR+5*chm.range_gate,1,'last')),profile(1:find(chm.range<=maxR+5*chm.range_gate,1,'last')),11,2);
   
   grad_profile_u = -grad_profile;grad_profile_u = grad_profile_u + nanmax(abs(grad_profile_u));grad_profile_u = grad_profile_u(~isnan(grad_profile_u));
   u = invariant_probability_distribution_silagadze(grad_profile_u,5);u = [NaN(5,1);u;NaN(5,1)];
   
%    u = -grad_profile;

   [~,ind_max] = local_maxima(u);
   val_max = grad_profile(ind_max);
   ind_max = ind_max(val_max<0);
   val_max = val_max(val_max<0);
   if isempty(val_max)
       continue;
   end
   val_max = val_max(chm.range(ind_max)>=minR & chm.range(ind_max)<=maxReffective & grad_profile(ind_max)<=log10(0.99)/(2*chm.range_gate));
   ind_max = ind_max(chm.range(ind_max)>=minR & chm.range(ind_max)<=maxReffective & grad_profile(ind_max)<=log10(0.99)/(2*chm.range_gate));
   if isempty(val_max)
       continue;
   end
   
   i_last = find(val_max<=log10(0.85)/(2*chm.range_gate),1,'first');
   if ~isempty(i_last)
      val_max = val_max(1:i_last);
      ind_max = ind_max(1:i_last);
   end
   
   cluster = {};
   cluster{end+1} = [ind_max(1),val_max(1)];
   for i=2:length(ind_max)
       if ind_max(i)-ind_max(i-1)<11
          cluster{end}(end+1,:) = [ind_max(i),val_max(i)];
       else
           cluster{end+1} = [ind_max(i),val_max(i)];
       end
   end
   ind_max = NaN(length(cluster),1);
   val_max = NaN(length(cluster),1); 
   for k=1:length(cluster)
       [~,indmincluster] = nanmax(cluster{k}(:,2));
       ind_max(k) = cluster{k}(indmincluster,1);
       val_max(k) = cluster{k}(indmincluster,2);
   end
   
   gradients(1,j+1) = chm.range(ind_max(1));
   
   [~,isort] = sort(val_max);
   r_sort = chm.range(ind_max(isort));
   if length(r_sort)>=3
       gradients(2,j+1) = r_sort(end-2);
   end
   if length(r_sort)>=2
       gradients(3,j+1) = r_sort(end-1);
   end
   gradients(4,j+1) = r_sort(end);
   
   if do_plot
       
       subplot(1,2,1)
       plot(profile,chm.range,'.-');
       xlim([3.5 6]);
       ylim([0 maxR]);
       set(gca,'ytick',0:100:maxR);
       grid on;box on;
       hold on;
       plot(xlim,minR*ones(1,2),'--k');
       plot(xlim,rSNR*ones(1,2),'--r');
       plot(xlim,nanmin(cbh)*ones(1,2),'--','color',0.5*ones(1,3));
       for l=1:length(ind_max)
           plot(profile(ind_max),chm.range(ind_max),'ok','markerfacecolor','k');
       end
       hold off;
       title(datestr(floor(chm.time(1))+xtime(j+1)));
       
       subplot(1,2,2)
       plot(grad_profile,chm.range(1:find(chm.range<=maxR+5*chm.range_gate,1,'last')),'.-');
       xlim([log10(0.5)/(2*chm.range_gate) -log10(0.5)/(2*chm.range_gate)]);
       ylim([0 maxR]);
       set(gca,'ytick',0:100:maxR);
       grid on;box on;
       hold on;
       plot(xlim,minR*ones(1,2),'--k');
       plot(xlim,rSNR*ones(1,2),'--r');
       plot(xlim,nanmin(cbh)*ones(1,2),'--','color',0.5*ones(1,3));
       for l=1:length(ind_max)
           plot(grad_profile(ind_max),chm.range(ind_max),'ok','markerfacecolor','k');
       end
       hold off;
       title(datestr(floor(chm.time(1))+xtime(j+1)));
       pause;
       
   end
   
end
end