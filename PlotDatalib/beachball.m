function handle = beachball(strike,dip,rake,x0,y0,radius,color,handle) 
% Usage: handle = beachball(strike,dip,rake,x0,y0,radius,color,handle) 
% 
% Plot a lower-hemisphere focal mechanism centered at x0,y0 
% with radius radius. 
% handle is an optional argument. If specified, and if handle 
% points to an existing beachball, then the beachball is updated. 
% Otherwise a new beachball is created and its handle is returned. 
    
 
   nargs = nargin; 
   if nargs < 3 
      error('Not enough args for beachball!'); 
   end 
 
   if nargs == 3 
      x0 = 0; 
      y0 = 0; 
      radius = 1; 
      color = 'k'; 
      handle = []; 
   end 
    
   if nargs == 4 
      error('Must specify either both x0 and y0 or neither of them!'); 
   end 
    
   if nargs == 5 
      radius = 1; 
      color = 'k'; 
      handle = []; 
   end 
 
   if nargs == 6 
      color = 'k'; 
      handle = []; 
   end 
 
   if nargs == 7 
      handle = []; 
   end 
     
       
   if radius <= 0 
      error('Radius must be positive!') 
   end 
    
   if isempty(handle) 
      handle = CreateBeachBall(strike,dip,rake,x0,y0,radius,color); 
   else 
      ModifyBeachBall(strike,dip,rake,x0,y0,radius,color,handle); 
      ModifyBeachBall(strike,dip,rake,x0,y0,radius,color,handle); 
   end 
    
% ------------------------------------------------------------------------------ 
 
 
 
 
 
function ModifyBeachBall(strike,dip,rake,x0,y0,radius,color,handle) 
 
 
   [x1,y1,x2,y2, Xp, Yp] = Boundaries(strike,dip,rake,x0,y0,radius); 
   set(handle.patch1,'Xdata',x1,'Ydata',y1); 
   set(handle.patch1,'FaceColor',color); 
    
   set(handle.patch2,'Xdata',x2,'Ydata',y2); 
   set(handle.patch2,'FaceColor',color); 
    
   azimuth = (0:360) *pi / 180; 
   x = x0 + cos(azimuth) * radius; 
   y = y0 + sin(azimuth) * radius; 
   set(handle.Equator,'Xdata',x,'Ydata',y); 
    
   set(handle.Paxis1,'Position',[Xp(1), Yp(1)]); 
   set(handle.Paxis2,'Position',[Xp(2), Yp(2)]); 
    
% -------------------------------------------------------------------------- 
 
 
    
    
    
 
function handle = CreateBeachBall(strike,dip,rake,x0,y0,radius,color) 
% Draw focal mechanism plot for current strike, dip, rake  
 
 
   [x1,y1,x2,y2,Xp, Yp] = Boundaries(strike,dip,rake,x0,y0,radius); 
   handle.patch1 = patch(x1,y1,color,'erasemode','xor',... 
                   'Tag','Patch1'); 
   handle.patch2 = patch(x2,y2,color,'erasemode','xor',... 
                   'Tag','Patch2'); 
   azimuth = (0:360) *pi / 180; 
   x = x0 + cos(azimuth) * radius; 
   y = y0 + sin(azimuth) * radius; 
   hBdr = line(x,y,'color','k','erasemode','background', 'Tag','Equator'); 
    
   handle.Equator   = hBdr; 
 
   handle.Paxis1 = text(Xp(1), Yp(1), 'P', 'color','k','erasemode','background',... 
                  'Tag', 'Paxis1','HorizontalAlignment','center',...
                  'VerticalAlignment', 'middle','fontsize',8); 
   handle.Paxis2 = text(Xp(2), Yp(2), 'P', 'color','k','erasemode','background',... 
                  'Tag', 'Paxis2','HorizontalAlignment','center',... 
                  'VerticalAlignment', 'middle','fontsize',8);
% ---------------------------------------------------------------------------- 
 
 
 
 
 
function [x1,y1,x2,y2, Xp, Yp] = Boundaries(strike,dip,rake,x0,y0,radius) 
 
% Get the boundaries of the compressional quadrants by starting with a  
% normalized fault (strike = 0, dip = 90, rake = 0) 
% Rotate the 1st and third quadrants of this system to the actual fault 
% orientation and then project onto equal-area lower hemisphere. 
 
   R = rotationMatrix(strike,dip,rake); 
 
 
   conv = pi/180; 
   angle = (0:180) * conv; 
   angle = angle(:)'; 
 
   SI = sin(angle); 
   ZE = zeros(size(angle)); 
   CS = cos(angle); 
 
% get projection of equatorial plane on normalized system 
   th2 = (0:360)*conv; 
   xb = cos(th2); 
   yb = sin(th2); 
   VV = [xb;yb;zeros(size(xb))]; 
   EqPlane = inv(R) * VV; 
 
 
 
% plane 1 
   V = [SI; ZE; CS];        %create 1/2 circle in +x-z plane 
   [xp1,yp1] = GetProjection(V, R); 
 
% plane 2 
   V = [ZE; SI; CS];        %create 1/2 circle in y-z plane 
   [xp2,yp2] = GetProjection(V, R); 
 
 
%  compressional part of equatorial plane connecting plane1 and plane2 
   II = find(EqPlane(1,:) >=0 &EqPlane(2,:) >=0); 
   VV=EqPlane(:,II); 
 
   [xxe,yye] = GetProjection2(VV,R); 
 
   [xp,yp] = Join(xp1,yp1,xp2,yp2,xxe, yye); 
   x1 = radius * xp + x0; 
   y1 = radius * yp + y0; 
 
 
 
    
 
% plane 3 
   V = [-SI; ZE; CS];        %create 1/2 circle in -x-z plane 
   [xp3,yp3] = GetProjection(V, R); 
    
% plane 4 
   V = [ZE; -SI; CS];        %create 1/2 circle in -y-z plane 
   [xp4,yp4] = GetProjection(V, R); 
 
%  compressional part of equatorial plane connecting plane3 and plane4 
 
   II = find(EqPlane(1,:) <=0 &EqPlane(2,:) <=0); 
   VV=EqPlane(:,II); 
   [xxe,yye] = GetProjection2(VV,R); 
   [xxp,yxp] = Join(xp3,yp3,xp4,yp4,xxe,yye); 
    
   x2 = radius * xxp + x0; 
   y2 = radius * yxp + y0; 
 
% Get projection of P-axis 
   Paxis = [-1 1;1 -1;0 0] /sqrt(2); 
   [Xpaxis, Ypaxis] = GetProjection(Paxis, R); 
   Xp = Xpaxis * radius + x0; 
   Yp = Ypaxis * radius + y0; 
    
% This must always be a 2-element vector even when only one pole is displayed    
   if length(Xp) == 1 
      Xp(2) = 1000; 
      Yp(2) = 1000; 
   end 
       
% ----------------------------------------------------------------       
 
 
 
 
 
function [xp,yp] = Join(xp1,yp1,xp2,yp2,eqx,eqy) 
   xp = []; 
   yp = []; 
   N = length(xp1); 
   M = length(xp2); 
   L = length(eqx); 
 
   % First join the two fault planes forcing the joint at the 
   % endpoints of smallest radius 
   r = sqrt(xp1.^2 + yp1.^2); 
   if r(1) > r(N) 
      xp = xp1(:); 
      yp = yp1(:); 
   else 
      xp = flipud(xp1(:)); 
      yp = flipud(yp1(:)); 
   end 
 
 
   r = sqrt(xp2.^2 + yp2.^2); 
   if r(1) > r(M) 
      xp = [xp; flipud(xp2(:))]; 
      yp = [yp; flipud(yp2(:))]; 
   else 
      xp = [xp; xp2(:)]; 
      yp = [yp; yp2(:)]; 
   end 
    
   if isempty(eqx) 
      return 
   end 
 
% sometimes eqx-eqy comes in as a closed curve, so check endpoints and 
% remove last if necessary 
 
   az = atan2(eqy,eqx); 
   II1 = find(az >=0 & az < pi/2); 
   II2 = find(az >= pi/2 & az < pi); 
   II3 = find(az < -pi/2 & az >= -pi); 
   II4 = find(az < 0 & az >= -pi/2); 
   if isempty(II1) | isempty(II4) 
      az(II3) = 2*pi + az(II3); 
      az(II4) = 2*pi + az(II4); 
   end 
   [az,II] = sort(az); 
   eqx = cos(az); 
   eqy = sin(az); 
    
   r = sqrt( (eqx - xp(1)).^2 + (eqy - yp(1)).^2); 
   if r(1) > r(L) 
      xp = [xp; eqx(:)]; 
      yp = [yp; eqy(:)]; 
   else 
      xp = [xp; flipud(eqx(:))]; 
      yp = [yp; flipud(eqy(:))]; 
   end 
    
% ---------------------------------------------------------------    
 
 
 
function [xp,yp] = GetProjection(V, R) 
 
   xp = []; 
   yp = [];    
   VP = R * V;              %rotate to strike-dip-rake 
   I = find(VP(3,:) >= 0);  %select part of rotated plane with + z 
   VPP = VP(:,I); 
   if isempty(VPP),return;end 
    
   r = sqrt(VPP(1,:).^2 + VPP(2,:).^2); 
   inc = ones(size(r)) * pi/2; 
   II = find(VPP(3,:) ~= 0); 
   if ~isempty(II) 
      inc(II) = atan(r(II) ./ VPP(3,II) ); 
   end 
 
   thet = atan2(VPP(2,:) , VPP(1,:)); 
   R0 = sqrt(2) * sin(inc/2); 
   xp = R0 .* sin(thet); 
   yp = R0 .* cos(thet); 
% ---------------------------------------------------------------- 
 
 
 
 
 
 
 
 
 
function [xp,yp] = GetProjection2(V, R) 
   % These points are guaranteed to be on the equator... 
   xp = []; 
   yp = [];    
   VP = R * V;              %rotate to strike-dip-rake 
   if isempty(VP),return;end 
 
   thet = atan2(VP(2,:) , VP(1,:)); 
   R0 = 1; 
   xp = R0 .* sin(thet); 
   yp = R0 .* cos(thet); 
% ---------------------------------------------------------------- 
 
 
 
 
 
 
 
    
function R = rotationMatrix(strike,dip,rake) 
 
   conv = pi/180; 
   phi = strike * conv; 
   delta = -(90 - dip) * conv; 
   lambda = rake * conv; 
   cp = cos(phi); 
   sp = sin(phi); 
   cd = cos(delta); 
   sd = sin(delta); 
   cl = cos(lambda); 
   sl = sin(lambda); 
   R3 = [cp -sp 0;sp cp 0; 0 0 1];   % rotation around Z for strike 
   R2 = [1 0 0 ; 0 cd -sd; 0 sd cd]; % rotation around X for dip 
   R1 = [cl 0 sl; 0 1 0; -sl 0 cl];  % rotation around Y for rake 
   R  = R3*R2*R1; 
% -----------------------------------------------------------------        
 
 
 
 
