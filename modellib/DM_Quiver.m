function [H,title]=DM_Quiver(XY,d,dcov,scale,origin)
%
%Adds velocity vector and error ellipses to the map

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Sep 9, 2000   Peter                 Major re-write
%   April, 2000   Peter Cervelli		    Original Code
%
%-------------------------------------------------------------

ns=size(XY,2);

%Arrow shafts

   switch length(d)==ns

      case 0 %Horizontal

         U=d(1:3:end)'*scale*1e2;
         V=d(2:3:end)'*scale*1e2;
         title='Horizontal GPS velocity and error ellipses';

      case 1 %Vertical

         V=d'*scale*1e2;
         U=zeros(size(V));
         title='Vertical GPS velocity and error ellipses';
   end

   X=[XY(1,:);XY(1,:)+U;repmat(NaN,1,ns)];
   Y=[XY(2,:);XY(2,:)+V;repmat(NaN,1,ns)];

%Arrow heads

   alpha=0.2;
   beta=0.33;
   Up=U;
   Vp=V;
   L=sqrt(sum([X(1,:)-X(2,:);Y(1,:)-Y(2,:)].^2));
   I=find(L>3);
   Up(I)=Up(I)./(L(I)/3);
   Vp(I)=Vp(I)./(L(I)/3);
   X=[X;X(2,:)-alpha*(Up+beta*(Vp+eps));X(2,:);X(2,:)-alpha*(Up-beta*(Vp+eps));repmat(NaN,1,ns)];
   Y=[Y;Y(2,:)-alpha*(Vp-beta*(Up+eps));Y(2,:);Y(2,:)-alpha*(Vp+beta*(Up+eps));repmat(NaN,1,ns)];

%Error ellipses

   if ~isempty(dcov)
      switch length(d)==ns

         case 0 %Horizontal

            r=linspace(0,2*pi,100);
            for i=ns:-1:1
               I=(i-1)*3+1;
               DCOV=scale^2*1e4*full(dcov(I:I+1,I:I+1));
               [v,w]=eig(DCOV);
               az=pi/2-atan2(v(1,1),v(2,1));
               w=sqrt(diag(w)*5.9915);
               elpts=[cos(az) -sin(az);sin(az) cos(az)]*[w(1)*cos(r);w(2)*sin(r)];
               X(8:108,i)=[elpts(1,:)'+X(2,i);NaN];
               Y(8:108,i)=[elpts(2,:)'+Y(2,i);NaN];
            end
            
         case 1 %Vertical

      end

   end

%Change to LLH if necessary

   if nargin ==5
      LLH=DM_local2llh([X(:)';Y(:)'],origin);
      X=LLH(1,:);
      Y=LLH(2,:);
   end

%Add vector field to map

   H=line('XData',        X(:), ...
          'YData',        Y(:), ...
          'ZData',        zeros(prod(size(X)),1), ...
          'Tag',          'DataVector', ...
          'LineWidth',    1, ...
          'HitTest',      'off');

