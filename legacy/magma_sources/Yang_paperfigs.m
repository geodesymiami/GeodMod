clear all ; close all

nu=0.25    % Poisson's ratio

f=figure('position',[259 75 960 768],'paperposition',[0.5 1.0 25.0 21.0],'inverthardcopy','off','Color',[1 1 1]);

% Figure 4a and 4b:  Influence of plunge variations
      y1=[-4:0.1:4];
      x1=zeros(1,length(y1));
	  %q=x1; x1=y1; y1=q;

      param_yang =[0; 0; 1.0 ;  1.0 ; 0.5 ; 0.25 ; 0.01 ; 0.01 ] ; % x y depth strength major-ax minor-ax strike plunge
      [u_yang1]=yang_source([param_yang(1:7);  0.01] , [x1 ; y1],nu); 
      [u_yang2]=yang_source([param_yang(1:7); 45.01] , [x1 ; y1],nu); 
      [u_yang3]=yang_source([param_yang(1:7); 89.08] , [x1 ; y1],nu); 

      axes
      plot(y1,u_yang1(3,:),'k-',y1,u_yang2(3,:),'r-',y1,u_yang3(3,:),'b-') ; axis tight
	  legend('th=0','th=45','th=90'); title('Fig. 4a U_v  z=1  a=0.5 c=0.25');
      set(gca,'NextPlot','add','position',[0.1 0.70 0.26 0.25])

      axes
      plot(y1,u_yang1(2,:),'k-',y1,u_yang2(2,:),'r-',y1,u_yang3(2,:),'b-') ; axis tight
	  legend('th=0','th=45','th=90',4); title('Fig. 4b  U_h');
      set(gca,'NextPlot','add','position',[0.4 0.70 0.26 0.25])

% Figure 5a and 5b:  Influence of depth (plunge=90)
      y1=[-2:0.1:2]; x1=zeros(1,length(y1));
      param_yang = [0; 0; 1.0 ;  1.0 ; 0.45 ; 0.5 ; 0.01 ; 89.99 ] ; % x y depth strength major-ax minor-ax strike plunge
      param_yang(5)=0.4 ; param_yang(6)=0.9*param_yang(5) ; [u_yang1]=yang_source(param_yang, [x1 ; y1],nu); 
      param_yang(5)=0.6 ; param_yang(6)=0.9*param_yang(5) ; [u_yang2]=yang_source(param_yang, [x1 ; y1],nu); 
      param_yang(5)=0.8 ; param_yang(6)=0.9*param_yang(5) ; [u_yang3]=yang_source(param_yang, [x1 ; y1],nu); 
      axes
      plot(y1,u_yang1(3,:),'k-',y1,u_yang2(3,:),'r-',y1,u_yang3(3,:),'b-') ; axis tight
	  legend('a/z=0.4','a/z=0.6','a/z=0.8'); title('Fig. 5a  U_v   z=1,c=0.9a,theta=90 ');
      set(gca,'NextPlot','add','position',[0.1 0.38 0.26 0.25])
      axes
      plot(y1,u_yang1(2,:),'k-',y1,u_yang2(2,:),'r-',y1,u_yang3(2,:),'b-') ; axis tight
	  legend('a/z=0.4','a/z=0.6','a/z=0.8',4); title('Fig. 5a  U_h   z=1,c=0.9a,theta=90 ');
      set(gca,'NextPlot','add','position',[0.4 0.38 0.26 0.25])

param_yang =[0; 0; 1.0 ;  1.0 ; 0.45 ; 0.5 ; 0.01 ; 0.01 ] ; % x y depth strength major-ax minor-ax strike plunge
% Figure 8a:

	  x1=[0.0:0.2:10] ; y1=zeros(1,length(x1)); 
      param_yang =[0; 0; 2.0 ;  1.0 ; 0.80 ; 0.28 ; 0.01 ; 89.9 ] ; % x y depth strength major-ax minor-ax strike plunge
      [u_yang]=yang_source(param_yang , [x1 ; y1],nu);
	  u_yang=u_yang./max(u_yang(3,:)) ;
      axes
      x=x1./param_yang(3);
      plot(x,u_yang(3,:),'k--',x,u_yang(1,:),'r--',x,u_yang(3,:),'k',x,u_yang(1,:),'r') ; axis tight
	  legend('Uv','Ur'); title('Fig. 8a');
      set(gca,'NextPlot','add','position',[0.10 0.03 0.26 0.25])

% Figure 8b:

	  x1=[0.0:0.2:5] ; y1=zeros(1,length(x1)); 
      param_yang =[0; 0; 1.0 ;  1.0 ; 0.80 ; 0.28 ; 0.01 ; 89.9 ] ; % x y depth strength major-ax minor-ax strike plunge
      [u_yang]=yang_source(param_yang , [x1 ; y1],nu);
	  u_yang=u_yang./max(u_yang(3,:)) ;
      axes
      x=x1./param_yang(3);
      plot(x,u_yang(3,:),'k--',x,u_yang(1,:),'r--',x,u_yang(3,:),'k',x,u_yang(1,:),'r') ; axis tight
	  legend('Uv','Ur'); title('Fig. 8b');

      set(gca,'NextPlot','add','position',[0.4 0.03 0.26 0.25])

