function  [OutVel,Parameters] = lockedandcreep(fault_para,dap)


%
% Function to model creep and locked depth. It generates displacement and
% strain for a screw dislocation.
%
% OutVel = lockedandcreep(fault_para DistAlongProf)
%
%  Inputs: (can be multidimentionnal (output of ngrid))
%
%  - fault_para = [ld;dsd;usd;loff;ffv] as:
%      ld   : Locking Depth                        (km) 
%      dsd  : Downward Slipping Depth              (km)
%      usd  : Upward Slipping Depth                (km) 
%      loff : Position offsets [offmin:inc:offmax] (km)
%      ffv  : Far Field Velocity                   (m/yr) 
%
%  - dap  : Distance along profile 	               (km)
%
%
%  Outputs:
%
%   OutVel     : Velocity
%   Parameters : Matrices with parameters combination
%

% N. Gourmelen modified from Unknown source April 2007
% Modified 01/2008 N. Gourmelen: added matrix computation and offset in fault position 	
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start modelling %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ld = fault_para(1,:);  dsd = fault_para(2,:);  usd = fault_para(3,:); loff = fault_para(4,:); ffv = fault_para(5,:);  
% conversion from km to meter of location parameter
  ld   = ld   * 1000;
  dsd  = dsd  * 1000;
  usd  = usd  * 1000;
  loff = loff * 1000;
  dap  = dap  * 1000;

  dap = dap(:);  dsd(find(dsd==0)) = 1e-10;  usd(find(usd==0)) = 1e-11;

% Combine inputs if multidimensionnals

if sum(size(ld)+size(ffv)+size(dsd)+size(usd)) > 8
	%[ldn,ffvn,dsdn,usdn]=ndgrid(ld,ffv,dsd,usd);
    exclude = find(dsd <= usd);  ld(exclude) = [];  ffv(exclude) = [];  dsd(exclude) = [];  usd(exclude) = []; % Lower depth for creeping should be higher than upper
    exclude = find(ld  <= dsd);  ld(exclude) = [];  ffv(exclude) = [];  dsd(exclude) = [];  usd(exclude) = []; % Locking depth should be higher that lower creeping depth
	ffv = ffv(:);  ld = ld(:);  dsd = dsd(:);  usd = usd(:);
    Parameters=[ffv(:) ld(:) dsd(:) usd(:)];
    %clear ffvn ldn dsdn usdn exclude
end

dap_saved = dap;

for ni=1:length(loff)

	dap = dap_saved+loff(ni);

	% Shallow creep between depths dsd and usd

	v1a = (ffv/pi);  v1an = repmat(v1a,1,length(dap));  
	v1b = ( atan( repmat(dap,1,length(usd)) ./ repmat(usd,1,length(dap))' ) - atan( repmat(dap,1,length(dsd)) ./ repmat(dsd,1,length(dap))' ) )';
	v1  = v1b .* v1an;

	% Free slip below ld

	v2a = (ffv/pi);  v2an = repmat(v2a,1,length(dap));  v2b = atan(repmat(dap,1,length(ld)) ./ repmat(ld,1,length(dap))')';
	v2  = v2b .* v2an;

	% Combine

	OutVel(:,:,ni) = (v1+v2)';
end
