function [igram]=add_mask_to_igram(igram,mask,lopt,opt)
%add_mask_to_igram      - multiplies igram.data with a mask
% usage:   [igram]=add_mask_to_igram(igram,mask,lopt)
%
%  Part of the TimeSeries suite
%  FA, March 2005,

if isfield(lopt,'subset')
   if size(mask)~=size(igram(1).data)        %subset mask only if igram and mask have different sizes
      mask=mask(lopt.subset(2):lopt.subset(2)+lopt.subset(4)-1,lopt.subset(1):lopt.subset(1)+lopt.subset(3)-1);    
   end
end

for i=1:length(igram)
    igram(i).data=igram(i).data.*mask;
end

% Add a phase ambiguity to certain area. Requires a dataset of same size as igram
% which contains values > 0 in are 
% can be created using the extractregion.m tool, 
if(opt.ambiguityOffset~=0)
           S=load(opt.ambiguityOffsetFile);
          offsetmask = S.offsetmask;
          offsetmask(find(isnan(offsetmask)))  = NaN;
          offsetmask(find(offsetmask==0))      = 0;
          offsetmask(find(offsetmask>0))       = 6.2832 * opt.ambiguityOffset;% 2pi
            for i=1:length(igram)
            igram(i).data=igram(i).data + offsetmask;
            end
          
end
%%% Tini
return
