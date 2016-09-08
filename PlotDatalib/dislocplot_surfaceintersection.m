function  dislocplot_surfaceintersection(disgeom,x_unit,varargin) 
%   displot        - Plot surface projection of the dislocation.
%    Input: Dislocation Geometry in a Global E-N-Up coordinate system
%		disgeom(1) = len = fault length in strike direction (km)
%		disgeom(2) = wid = fault width in dip direction (km)
%		disgeom(3) = dep = fault depth (km)
%		disgeom(4) = dip = dip angle (degrees)
%		disgeom(5) = strik =  strike, clockwise from N (degrees)
%		disgeom(6) = delE  = East offset of midpoint from origin (km)
%		disgeom(7) = delN  = North offset of midpoint from origin (km)
%		disgeom(8) = ss  =  strike slip motion (m)	
%		disgeom(9) = ds  =  dip slip motion (m)	
%		disgeom(10) = op  =  opening  motion (m)	
%
%		x_unit : km (defaults) or (degrees,lola)	
%
%		May 2007, based on Paul Segall's displot which works only in km
%       7/2008 FA  For ll changed right,left so that compatible with
%       right,left for xy

 switch length(varargin)
 case 1
        PlotType = varargin{1};
 case 2
        Profile  = varargin{1};
        PlotType = varargin{2};
 end

 i = sqrt(-1);
%
if exist('x_unit','var') && (strcmp(x_unit,'degres') | strcmp(x_unit,'degrees'))
  
  len          = disgeom(1); 
  wid          = disgeom(2); 
  dep          = disgeom(3); 
  delE         = disgeom(6); 
  delN         = disgeom(7);
  strike       = disgeom(5);
  strkr        = (90-disgeom(5))*pi/180;
  dip          = disgeom(4);
  dipr         = disgeom(4)*pi/180;
    
  Horiz_shift  = dep / tan( dip/180*pi );
  
  uc_surfaceintersect = reckon(disgeom(7),disgeom(6),km2deg(Horiz_shift), strike+270) ; % center of surface interection
  
  % center of upper edge
  %uc =  reckon(disgeom(7),disgeom(6),km2deg( wid*cos(dipr)), strike+270) ;
  %
  uc = uc_surfaceintersect;
  lr =  reckon(disgeom(7),disgeom(6),km2deg( 0.5*len), strike) ;
  ll =  reckon(disgeom(7),disgeom(6),km2deg(-0.5*len), strike) ;
  ur =  reckon(uc(1)     ,uc(2)     ,km2deg( 0.5*len), strike) ;
  ul =  reckon(uc(1)     ,uc(2)     ,km2deg(-0.5*len), strike) ;

  %plot([ul(2),ll(2),lr(2),ur(2)],[ul(1),ll(1),lr(1),ur(1)]);% , axis('equal');
  plot([ul(2),ur(2)],[ul(1),ur(1)],'g')                    ;% , axis('equal');
else
  
  disgeom = Extend2Surface(disgeom);
  len     = disgeom(1); 
  wid     = disgeom(2); 
  dep     = disgeom(3); 
  delE    = disgeom(6); 
  delN    = disgeom(7);
  strike  = disgeom(5);
  strkr   =(90-disgeom(5))*pi/180;
  dip     = disgeom(4);
  dipr    = disgeom(4)*pi/180;

  ll      = -0.5*len +i*0;  lr = 0.5*len + i*0; 
  ul      = -0.5*len +i*wid*cos(dipr); ur = 0.5*len + i*wid*cos(dipr);

  %transform corners to global coordinates
  strkr   = (90-disgeom(5))*pi/180;

  ll      = (delE+delN*i) + ll*exp(i*strkr);
  lr      = (delE+delN*i) + lr*exp(i*strkr);
  ul      = (delE+delN*i) + ul*exp(i*strkr);
  ur      = (delE+delN*i) + ur*exp(i*strkr);

  %plot( [real(ul), real(ur)], [imag(ul), imag(ur)], 'g' ), axis(axis);
  switch PlotType
  case {'2D'}
       plot( [real(ll), real(lr)], [imag(ll), imag(lr)],'g' ), axis(axis);
  case {'1D'}
       disloc_profile_dist = mean( Profile.unitvec*[ [real(ll) real(lr)]-Profile.xy(1,1); [imag(ll), imag(lr)]-Profile.xy(1,2) ] ); % mean because Profile may not be perpendicular to dilocation
       plot(disloc_profile_dist,-0.1,'^k','LineWidth',4);
  end
end
