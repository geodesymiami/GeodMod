function  displot(disgeom) 
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
%%
len =disgeom(1); wid = disgeom(2); dep = disgeom(3); delE=disgeom(6); delN=disgeom(7);
i = sqrt(-1);
dipr = disgeom(4)*pi/180;
%
% corners in local coordinates
ll = -0.5*len +i*0;  lr = 0.5*len + i*0; 
ul = -0.5*len +i*wid*cos(dipr); ur = 0.5*len + i*wid*cos(dipr);
%
%transform corners to global coordinates
strkr = (90-disgeom(5))*pi/180;
%
ll = (delE+delN*i) + ll*exp(i*strkr);
lr = (delE+delN*i) + lr*exp(i*strkr);
ul = (delE+delN*i) + ul*exp(i*strkr);
ur = (delE+delN*i) + ur*exp(i*strkr);
%
plot( [real(ll), real(lr)], [imag(ll), imag(lr)],'LineWidth',2 ), axis('equal');
plot( [real(ll), real(ul)], [imag(ll), imag(ul)],'LineWidth',2 ), axis(axis);
plot( [real(ul), real(ur)], [imag(ul), imag(ur)], 'g','LineWidth',2 ), axis(axis);
plot( [real(ur), real(lr)], [imag(ur), imag(lr)],'LineWidth',2 ), axis(axis);
