% Sjonni's old code.
% Note that the RMS is defined differently then in geodmod. I am not sure which one to use  (FA 7/08)

pred =  K_weighted*s ./sqrt(normalization);

resi = dprime-pred;
resi = resi.*sqrt(normalization);

% need to get rid of the following commands %sigphi_1 = dataset(1).sigphi;sigphi_2= dataset(2).sigphi; sigphi_3 = dataset(3).sigphi;
readfrom_dataset_structure;      
%inverseopt.distribopt.slip             = [0 0 1]
%data_1  = dataset(1).Ndata ; Ndata_2 = dataset(2).Ndata ; Ndata_3 = dataset(3).Ndata ;
%DataSet_1=dataset(1).DataSet;DataSet_2=dataset(2).DataSet;DataSet_3=dataset(3).DataSet; 

    tot_rms  = norm(resi)                            / sqrt(length(resi));
    rms_1    = norm(resi(1:datind(1)))                 / sqrt(Ndata_1);
if D_2 rms_2 = norm(resi(datind(1)+1:datind(2))) / sqrt(Ndata_2); end
if D_3 rms_3 = norm(resi(datind(2)+1:datind(3))) / sqrt(Ndata_3); end
if D_4 rms_4 = norm(resi(datind(3)+1:datind(4))) / sqrt(Ndata_4); end
if D_5 rms_5 = norm(resi(datind(4)+1:datind(5))) / sqrt(Ndata_5); end

        St = sprintf('RMS              : %2.2f cm',tot_rms*100); disp(St)
        St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_1,rms_1*100); disp(St)
if D_2  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_2,rms_2*100); disp(St) ; end
if D_3  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_3,rms_3*100); disp(St) ; end
if D_4  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_4,rms_4*100); disp(St) ; end
if D_5  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_5,rms_5*100); disp(St) ; end

% Calculate the weight of each data set
%Wper = [sum(W(1:Ndata_1)) sum(W(Ndata_1+1:Ndata_1+Ndata_2))]';
        W    = 1./sqrt(normalization);
        Wper = [       sum(W(          1:datind(1)))];
if D_2  Wper = [Wper ; sum(W(datind(1)+1:datind(2)))];  end ;
if D_3  Wper = [Wper ; sum(W(datind(2)+1:datind(3)))];  end ;
if D_4  Wper = [Wper ; sum(W(datind(3)+1:datind(4)))];  end ;
if D_5  Wper = [Wper ; sum(W(datind(4)+1:datind(5)))];  end ;
        Wper = Wper/sum(Wper)*100;
disp(' ')
       St = sprintf('Weight dataset %s: %2.1f%s',DataSet_1,Wper(1),'%'); disp(St)
if D_2 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_2,Wper(2),'%'); disp(St); end
if D_3 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_3,Wper(3),'%'); disp(St); end
if D_4 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_4,Wper(4),'%'); disp(St); end
if D_5 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_5,Wper(5),'%'); disp(St); end
