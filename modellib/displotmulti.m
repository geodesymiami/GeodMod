function [] = displotmulti(param,objfunc,modelopt,x_unit)
%   displotmulti    - plots one or more dislocations using displot
% usage:  [] = displotmulti(param);
%
% FA, last modified 30 Nov 2003 
% Yunjun, 2015-12-05: add circle()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotType = '2D';

switch  objfunc
case {'GenericObjectiveFunction','ManyDislocOneInflation'}

     N_disloc       = modelopt.N_disloc;
     N_mogi         = modelopt.N_mogi;
     N_penny        = modelopt.N_penny;
     N_mctigue      = modelopt.N_mctigue;
     N_pCDM         = modelopt.N_pCDM;
     N_yang         = modelopt.N_yang;
     N_multidisloc  = modelopt.N_multidisloc;
     N_visco1d      = modelopt.N_visco1d;
     N_squaredisloc = modelopt.N_squaredisloc;
     
     if N_multidisloc  
         param    = [param(1:N_disloc*10); multidislocpar2dislocpar(param(N_disloc*10+1:end),modelopt.multidislocopt,x_unit)];
         N_disloc = modelopt.multidislocopt.N_disloc + N_disloc;
     end

%      if N_disloc >= 1   plot(param(6),param(7),'k*'); displot(param(1:10)); dislocplot_surfaceintersection(param(1:10),x_unit,PlotType); axis tight; param(1:10)=[]; end
%      if N_disloc >= 2   plot(param(6),param(7),'k*'); displot(param(1:10)); dislocplot_surfaceintersection(param(1:10),x_unit,PlotType); axis tight; param(1:10)=[]; end
     if N_disloc >= 1;
         plot(param(6),param(7),'k*');
         displot(param(1:10));
         axis tight; param(1:10)=[];
     end     % Surface intersection should not be plotted by default, need to be build an option for that, Falk, 2015-11-04
     if N_disloc >= 2;   plot(param(6),param(7),'k*'); displot(param(1:10)); axis tight; param(1:10)=[]; end
     if N_disloc >= 3;   plot(param(6),param(7),'k*'); displot(param(1:10)); axis tight; param(1:10)=[]; end
     if N_disloc >= 4;   plot(param(6),param(7),'k*'); displot(param(1:10)); axis tight; param(1:10)=[]; end
     if N_disloc >= 5;
                       for i=5:N_disloc
                           if i<=15 || i>N_disloc-15 %display only 45 dislocations to save time
                              plot(param(6),param(7),'k*');  displot(param(1:10)); axis tight;  
                           end
                           param(1:10)=[];
                       end
     end                   
     %if N_squaredisloc >= 1;   plot(param(10),param(11),'k*'); displot([param(1);param(1:9)]); axis tight; param(1:9)=[]; end
     if N_squaredisloc >= 1;   plot(param(5),param(6),'k*','MarkerSize',5); displot([param(1);param(1:9)]); axis tight; param(1:9)=[]; end
     
     if N_mogi   >= 1;      plot(param(1),param(2),'k+','MarkerSize',9,'LineWidth',3);       param(1:4) =[]; end
     if N_mogi   >= 2;      plot(param(1),param(2),'k+');                                    param(1:4) =[]; end
     if N_mogi   >= 3;      plot(param(1),param(2),'k+');                                    param(1:5) =[]; end
     if N_penny  == 1;
                            if param(4)>3;  plot(param(1),param(2),'k+','MarkerSize',9,'LineWidth',3);  end
                            circle(param(1),param(2),param(4));                              param(1:5) =[];
     if N_pCDM   >= 1;      plot(param(1),param(2),'k+','MarkerSize',9,'LineWidth',3);       param(1:4) =[]; end
     end
     if N_mctigue== 1;      
                            if param(3)>3;  plot(param(1),param(2),'k+','MarkerSize',9,'LineWidth',3);  end
                            circle(param(1),param(2),param(4));                              param(1:5) =[];
     end
     if N_yang   == 1;      plot(param(1),param(2),'k+');                                    param(1:5) =[]; end
     if N_visco1d>= 1;      plot(param(6),param(7),'k*');  displot(param(1:10)); axis tight; param(1:10)=[]; end
end
end

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)

ang=0:0.005:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'k','LineWidth',2);

end


