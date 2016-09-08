switch PhaseRamp
case{'Ramp',1}
   dummy_plane=[-10 ; -0.1 ; -0.1 ];                            %add dummy plane because of nnls                                                        
   if InverseSign  dummy_plane=-1*dummy_plane;  end                                                                                                     
            plane_1=[ ones(Ndata_1,1)  coord_1(1,:)'  coord_1(2,:)']; datavec_1= datavec_1-(plane_1*dummy_plane)' ;                             
   if D_2   plane_2=[ ones(Ndata_2,1)  coord_2(1,:)'  coord_2(2,:)']; datavec_2= datavec_2-(plane_2*dummy_plane)' ; end                         
   if D_3   plane_3=[ ones(Ndata_3,1)  coord_3(1,:)'  coord_3(2,:)']; datavec_3= datavec_3-(plane_3*dummy_plane)' ; end                         
   if D_4   plane_4=[ ones(Ndata_4,1)  coord_4(1,:)'  coord_4(2,:)']; datavec_4= datavec_4-(plane_4*dummy_plane)' ; end                         
   if D_5   plane_5=[ ones(Ndata_5,1)  coord_5(1,:)'  coord_5(2,:)']; datavec_5= datavec_5-(plane_5*dummy_plane)' ; end                         
case{'Const',1}
   dummy_plane=[-10 ;  0.0 ;  0.0 ];         %add dummy plane because of nnls                                                                           
   if InverseSign  dummy_plane=-1*dummy_plane;  end                                                                                                     
            plane_1=[ ones(Ndata_1,1)  coord_1(1,:)'  coord_1(2,:)']; datavec_1= datavec_1-(plane_1*dummy_plane)' ;                             
   if D_2   plane_2=[ ones(Ndata_2,1)  coord_2(1,:)'  coord_2(2,:)']; datavec_2= datavec_2-(plane_2*dummy_plane)' ; end                         
   if D_3   plane_3=[ ones(Ndata_3,1)  coord_3(1,:)'  coord_3(2,:)']; datavec_3= datavec_3-(plane_3*dummy_plane)' ; end                         
   if D_4   plane_4=[ ones(Ndata_4,1)  coord_4(1,:)'  coord_4(2,:)']; datavec_4= datavec_4-(plane_4*dummy_plane)' ; end                         
   if D_5   plane_5=[ ones(Ndata_5,1)  coord_5(1,:)'  coord_5(2,:)']; datavec_5= datavec_5-(plane_5*dummy_plane)' ; end                         
case{false,1}
            plane_1=zeros(Ndata_1,3) ;                                                                                                                  
   if D_2   plane_2=zeros(Ndata_2,3) ;   end                                                                                                            
   if D_3   plane_3=zeros(Ndata_3,3) ;   end                                                                                                            
   if D_4   plane_4=zeros(Ndata_4,3) ;   end                                                                                                            
   if D_5   plane_5=zeros(Ndata_5,3) ;   end                                                                                                            
end  
if D_2  Gadd = [plane_1           zeros(Ndata_1,3);
                zeros(Ndata_2,3)  plane_2];          end ;
if D_3  Gadd = [plane_1           zeros(Ndata_1,3)     zeros(Ndata_1,3) ;
                zeros(Ndata_2,3)  plane_2              zeros(Ndata_2,3) ;
                zeros(Ndata_3,3)  zeros(Ndata_3,3)     plane_3        ] ;   end
if D_4  Gadd = [plane_1           zeros(Ndata_1,3)     zeros(Ndata_1,3)     zeros(Ndata_1,3) ;
                zeros(Ndata_2,3)  plane_2              zeros(Ndata_2,3)     zeros(Ndata_2,3) ;
                zeros(Ndata_3,3)  zeros(Ndata_3,3)     plane_3              zeros(Ndata_3,3) ;
                zeros(Ndata_4,3)  zeros(Ndata_4,3)     zeros(Ndata_4,3)     plane_4         ]; end ;
if D_5  Gadd = [plane_1           zeros(Ndata_1,3)     zeros(Ndata_1,3)     zeros(Ndata_1,3)  zeros(Ndata_1,3) ;
                zeros(Ndata_2,3)  plane_2              zeros(Ndata_2,3)     zeros(Ndata_2,3)  zeros(Ndata_2,3) ;
                zeros(Ndata_3,3)  zeros(Ndata_3,3)     plane_3              zeros(Ndata_3,3)  zeros(Ndata_3,3) ;
                zeros(Ndata_4,3)  zeros(Ndata_4,3)     zeros(Ndata_4,3)     plane_4           zeros(Ndata_4,3) ;
                zeros(Ndata_5,3)  zeros(Ndata_5,3)     zeros(Ndata_5,3)     zeros(Ndata_5,3)  plane_5          ]; end ;
