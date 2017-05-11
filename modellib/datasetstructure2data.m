function [d,coord,normalization,radarlook,datind,hgt,G_phaseramp,D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8,D_9]=datasetstructure2data(dataset,inverseopt);
%
% datasetstructure2data   -  converts dataset structure into data format for modelling       
%
% usage: function [d,coord,normalization,radarlook,datind,hgt,D_1,D_2,D_3,D_4,D_5]=datasetstructure2data(dataset);
%
% Input:   dataset           dislocation geometry
%
% Output:  d  ,coord
%

% shortcut if dataset(1) has already produced previously
if isfield(dataset(1),'fulldata')
   fulldata=dataset(1).fulldata ;
   [d,coord,normalization,radarlook,datind,hgt,G_phaseramp]=deal(fulldata.d,fulldata.coord,fulldata.normalization,fulldata.radarlook,fulldata.datind,fulldata.hgt,fulldata.G_phaseramp);
   [D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]       =deal(fulldata.D_1,fulldata.D_2,fulldata.D_3,fulldata.D_4,fulldata.D_5,fulldata.D_6,fulldata.D_7,fulldata.D_8);
else
   D_1=true;
   D_2=length(dataset)>=2 ;
   D_3=length(dataset)>=3 ;
   D_4=length(dataset)>=4 ;
   D_5=length(dataset)>=5 ;
   D_6=length(dataset)>=6 ;
   D_7=length(dataset)>=7 ;
   D_8=length(dataset)>=8 ;

   %%% FA 7/2008: below does the same thing but is much easier
   %tmp=cell(length(dataset),1);;
   %[tmp{:}]=deal(dataset.datavec);  d=cell2mat(tmp')';
   %[tmp{:}]=deal(dataset.coord);        coord=cell2mat(tmp');
   %[tmp{:}]=deal(dataset.radarlook);    radarlook=cell2mat(tmp'); 
   %[tmp{:}]=deal(dataset.Ndata);        Ndata =cell2mat(tmp);      datind=cumsum(Ndata);
   %[tmp{:}]=deal(dataset.sigphi);       sigphi=cell2mat(tmp);
   %if isfield(dataset,'hgt') [tmp{:}]=deal(dataset.hgt);          hgt=cell2mat(tmp')';  else hgt=[]; end
   
   d         = [dataset.datavec]';
   coord     = [dataset.coord];
   radarlook = [dataset.radarlook];
   Ndata     = [dataset.Ndata]';
   sigphi    = [dataset.sigphi]';
   datind    = cumsum(Ndata);
   
   if isfield(dataset,'hgt')  hgt = [dataset.hgt]';  else hgt=[]; end

   normalization=[];

   if D_1 NORM_1=isfield(dataset(1),'normalization') && ~ isempty(dataset(1).normalization); end
   if D_2 NORM_2=isfield(dataset(2),'normalization') && ~ isempty(dataset(2).normalization); end
   if D_3 NORM_3=isfield(dataset(3),'normalization') && ~ isempty(dataset(3).normalization); end
   if D_4 NORM_4=isfield(dataset(4),'normalization') && ~ isempty(dataset(4).normalization); end
   if D_5 NORM_5=isfield(dataset(5),'normalization') && ~ isempty(dataset(5).normalization); end
   if D_6 NORM_6=isfield(dataset(6),'normalization') && ~ isempty(dataset(6).normalization); end
   if D_7 NORM_7=isfield(dataset(7),'normalization') && ~ isempty(dataset(7).normalization); end
   if D_8 NORM_8=isfield(dataset(8),'normalization') && ~ isempty(dataset(8).normalization); end

   if D_1  if NORM_1 normalization=[normalization;dataset(1).normalization'*sigphi(1)]; else tmp = [sigphi(1).*ones(Ndata(1),1)]; normalization = tmp;               end; end; % FA 5/17: removed tmp = tmp.^2
   if D_2  if NORM_2 normalization=[normalization;dataset(2).normalization'*sigphi(2)]; else tmp = [sigphi(2).*ones(Ndata(2),1)]; normalization=[normalization;tmp] ;end; end; 
   if D_3  if NORM_3 normalization=[normalization;dataset(3).normalization'*sigphi(3)]; else tmp = [sigphi(3).*ones(Ndata(3),1)]; normalization=[normalization;tmp] ;end; end; 
   if D_4  if NORM_4 normalization=[normalization;dataset(4).normalization'*sigphi(4)]; else tmp = [sigphi(4).*ones(Ndata(4),1)]; normalization=[normalization;tmp] ;end; end; 
   if D_5  if NORM_5 normalization=[normalization;dataset(5).normalization'*sigphi(5)]; else tmp = [sigphi(5).*ones(Ndata(5),1)]; normalization=[normalization;tmp] ;end; end; 
   if D_6  if NORM_6 normalization=[normalization;dataset(6).normalization'*sigphi(6)]; else tmp = [sigphi(6).*ones(Ndata(6),1)]; normalization=[normalization;tmp] ;end; end; 
   if D_7  if NORM_7 normalization=[normalization;dataset(7).normalization'*sigphi(7)]; else tmp = [sigphi(7).*ones(Ndata(7),1)]; normalization=[normalization;tmp] ;end; end; 
   if D_8  if NORM_8 normalization=[normalization;dataset(8).normalization'*sigphi(5)]; else tmp = [sigphi(8).*ones(Ndata(8),1)]; normalization=[normalization;tmp] ;end; end; 
   [G_phaseramp]=MakeDesignMatrix(dataset,[],inverseopt);

end
