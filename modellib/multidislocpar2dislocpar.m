function [dpar]  = multidislocpar2dislocpar(mpar,multidislocopt,x_unit)
%multidislocpar2dislocpar  - returns modelvector for multiple dislocations
%based on multidisloc parameters
%
%ToDo: There may be problems with negative HorzOffsets and how c_HorzOffset is calculated in this case (in particular for 3 dislocations)

 i  = sqrt(-1);
 
 N_disloc          = multidislocopt.N_disloc;
 ind               = multidislocopt.ind;
 orig_indHorzOff   = multidislocopt.orig_indHorzOff;
 orig_indVertOff   = multidislocopt.orig_indVertOff;
 orig_parind       = multidislocopt.orig_parind;
 new_parind        = multidislocopt.new_parind;
 
 dpar              = repmat(mpar(1:10),N_disloc,1);
 dpar(new_parind)  = mpar(orig_parind);                        % here the multidisloc parameters are inserted

 [len1,wid1,dep1,dip1,strike1,xe1, xn1 ,ss1,ds1,op1]=deal(dpar(1), dpar(2), dpar(3), dpar(4), dpar(5), dpar(6), dpar(7), dpar(8), dpar(9), dpar(10));
 [len2,wid2,junk,dip2,strike2,junk,junk,ss2,ds2,op2]=deal(dpar(11),dpar(12),dpar(13),dpar(14),dpar(15),dpar(16),dpar(17),dpar(18),dpar(19),dpar(20));
 [ind_dep,ind_xy]  = deal(13,[16 17]);

 if N_disloc == 3
    [len3,wid3,junk,dip3,strike3,junk,junk,ss3,ds3,op3]=deal(dpar(21),dpar(22),dpar(23),dpar(24),dpar(25),dpar(26),dpar(27),dpar(28),dpar(29),dpar(30));
    [ind_dep,ind_xy]  = deal([13 23],[16 17 26 27]);
 end
 
  if isempty(orig_indVertOff) 
             VertOff = 0; 
  else
             VertOff = mpar(orig_indVertOff);
             if (VertOff  == -999)  
                 VertOff1 = wid1/2;
                 VertOff2 = wid2/2;
                 if N_disloc == 3  VertOff3=wid3/2; end
             else 
                 VertOff1 = VertOff(1)/2;
                 VertOff2 = VertOff(1)/2;
                 if N_disloc == 3  VertOff3=VertOff(2)/2; end
             end
  end

  if isempty(orig_indHorzOff)
             HorzOff = 0; 
  else
             HorzOff      = mpar(orig_indHorzOff);
             if (HorzOff  == -999)  
                 HorzOff1 = len1/2;
                 HorzOff2 = len2/2;
                 if N_disloc == 3  HorzOff3=len3/2; end
             else 
                 HorzOff1 = HorzOff(1)/2;                          % the index are confusing: the same offset applies for the first and second third dislocation (first field, HorzOff(1))
                 HorzOff2 = HorzOff(1)/2;                          % the index are confusing: HorzOff3 (for the third dislocation) is actually the second field in HorzOff, HorzOff(2)
                 if N_disloc == 3  HorzOff3=HorzOff(2)/2; end
             end
  end
  
 switch x_unit
     case {'degrees' 'degres'}
        % center of upper edge

        lc1    = [xn1 xe1];
       
        lr1    = reckon(lc1(1),lc1(2),km2deg( 0.5*len1), strike1) ;
        ll1    = reckon(lc1(1),lc1(2),km2deg(-0.5*len1), strike1) ;
        uc1    = reckon(lc1(1),lc1(2),km2deg( wid1*cosd(dip1)), strike1+270) ;   
        ur1    = reckon(uc1(1),uc1(2),km2deg( 0.5*len1),strike1) ;
        ul1    = reckon(uc1(1),uc1(2),km2deg(-0.5*len1),strike1) ;
        
        tmplr1 = reckon(   lc1(1),   lc1(2),km2deg(HorzOff1),strike1);                                    % first go from lc1  half dislocation length im strike1 direction
        tmp    = reckon(tmplr1(1),tmplr1(2),km2deg(HorzOff2),strike2);                                    % then  go another half dislocation length in stike2 direction
        tmplc2 = reckon(tmp(1),   tmp(2),   km2deg(VertOff1*cosd(dip1)+VertOff2*cosd(dip2)),strike1+270); % finally go horizontally according to the horizontal component of the along-dip offset
        
        lc2    = fliplr(tmplc2);    
        dep2   = dep1 - VertOff1*sind(dip1) - VertOff2*sind(dip2);
        
        multidisloc_dep = dep2;
        multidisloc_xy  = lc2;
        
        if N_disloc == 3      
                    tmplr2 = reckon(tmplc2(1),tmplc2(2),km2deg(HorzOff2),strike2);                                   % first go from lc2  half dislocation length im strike2 direction
                    tmp    = reckon(tmplr2(1),tmplr2(2),km2deg(HorzOff3),strike3);                                   % then  go another half dislocation length in stike2 direction 
                    tmplc3 = reckon(tmp(1),   tmp(2),   km2deg(VertOff2*cosd(dip2)+VertOff3*cosd(dip3)),strike3+270);% finally go horizontally according to the horizontal component of the along-dip offset
                    
                    lc3    = fliplr(tmplc3);    
                    dep3   = dep1 - VertOff1*sind(dip1) - VertOff2*sind(dip2) - VertOff3*sind(dip2);
        
                    multidisloc_dep = [ dep2 dep3 ];
                    multidisloc_xy  = [ lc2  lc3  ];
        end   
        
        %logmessage('+++PLOTTING NEED TO BE CHECKED+++');                
     case {'km'} 
         
        lc_orig = xe1 + i*xn1;
        
        ll1    = (-0.5*len1 + i*0)               * exp(i*deg2rad(90-strike1));  
        lr1    = ( 0.5*len1 + i*0)               * exp(i*deg2rad(90-strike1));
        ul1    = (-0.5*len1 + i*wid1*cosd(dip1)) * exp(i*deg2rad(90-strike1)); 
        ur1    = ( 0.5*len1 + i*wid1*cosd(dip1)) * exp(i*deg2rad(90-strike1));
        uc1    = (            i*wid1*cosd(dip1)) * exp(i*deg2rad(90-strike1)); 
 
        ll2    = (-0.5*len2 + i*0)               * exp(i*deg2rad(90-strike2));  
        lr2    = ( 0.5*len2 + i*0)               * exp(i*deg2rad(90-strike2));
        ul2    = (-0.5*len2 + i*wid2*cosd(dip2)) * exp(i*deg2rad(90-strike2)); 
        ur2    = ( 0.5*len2 + i*wid2*cosd(dip2)) * exp(i*deg2rad(90-strike2));
        
        tmpll1 = (-HorzOff1  + i*0)                * exp(i*deg2rad(90-strike1));  
        tmplr1 = ( HorzOff1  + i*0)                * exp(i*deg2rad(90-strike1));
        tmpul1 = (-HorzOff1  + i*wid1*cosd(dip1))  * exp(i*deg2rad(90-strike1)); 
        tmpur1 = ( HorzOff1  + i*wid1*cosd(dip1))  * exp(i*deg2rad(90-strike1));
 
        tmpll2 = (-HorzOff2  + i*0)               * exp(i*deg2rad(90-strike2));  
        tmplr2 = ( HorzOff2  + i*0)               * exp(i*deg2rad(90-strike2));
        tmpul2 = (-HorzOff2  + i*wid2*cosd(dip2)) * exp(i*deg2rad(90-strike2)); 
        tmpur2 = ( HorzOff2  + i*wid2*cosd(dip2)) * exp(i*deg2rad(90-strike2));

        c_VertOff  = i*(VertOff1*cosd(dip1) + VertOff2*cosd(dip2)) * exp(i*deg2rad(90-strike1));      % FA 3/2010: do we need strike2 here?
        c_HorzOff  = tmplr1 - tmpll2;
        
        lc2  = lc_orig + c_VertOff + c_HorzOff;        
        dep2 = dep1 - VertOff1*sind(dip1) - VertOff2*sind(dip2);
        lc2  = [real(lc2) imag(lc2)];
        
        multidisloc_dep = dep2;
        multidisloc_xy  = lc2;
        
        if N_disloc == 3
           ll3    = (-0.5*len3 + i*0)               * exp(i*deg2rad(90-strike3));  
           lr3    = ( 0.5*len3 + i*0)               * exp(i*deg2rad(90-strike3));
           ul3    = (-0.5*len3 + i*wid3*cosd(dip3)) * exp(i*deg2rad(90-strike3)); 
           ur3    = ( 0.5*len3 + i*wid3*cosd(dip3)) * exp(i*deg2rad(90-strike3));
           
           tmpll3 = (-HorzOff3  + i*0)               * exp(i*deg2rad(90-strike2));  
           tmplr3 = ( HorzOff3  + i*0)               * exp(i*deg2rad(90-strike2));
           tmpul3 = (-HorzOff3  + i*wid3*cosd(dip3)) * exp(i*deg2rad(90-strike3)); 
           tmpur3 = ( HorzOff3  + i*wid3*cosd(dip3)) * exp(i*deg2rad(90-strike3));
           
           c_VertOff  = c_VertOff + i* VertOff3*cosd(dip3) * exp(i*deg2rad(90-strike1));                  % For vertical offset
           c_HorzOff  = c_HorzOff + tmplr2 - tmpll3;          
        
           lc3  = lc_orig + c_VertOff + c_HorzOff;        
           dep3 = dep1 - VertOff3*sind(dip1);
           lc3  = [real(lc3) imag(lc3)];
        
           multidisloc_dep = [dep2 dep3];
           multidisloc_xy  = [lc2  lc3 ];
        end   

 end

 dpar(ind_dep)    = multidisloc_dep;
 dpar(ind_xy)     = multidisloc_xy;
 

 
