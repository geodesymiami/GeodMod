function [Focals]=MakeFocals(focalfile,focalfileformat,basemap)
%MakeFocals           - reads Focalsfile and saves s structure as matlab file
%
%  usage:  [Focals]=MakeFocals(focalfile,focalfileformat)
%
%  basemap is needed to switch from lola to xy
%

global dir_out

  if focalfile(1) CheckInOut(focalfile,''); end             % check whether quakefile exist

  focals_matfile = fullfile(dir_out,'focals.mat');

  if  CheckInOut('',focals_matfile) load(focals_matfile); return; end
%
% open and read quakefile
%

  fid = fopen(focalfile);
  switch focalfileformat
		 case{'Focal'}
			%C	  =	 textscan(fid, '%f %f %f %f %f %f %f %f','headerLines',0 );
            C	  =	 textscan(fid, '%f %f %f %f %f %f %f %f');
            C=cell2mat(C);
			lat	  =  [C(:,1)'];
			long  =  [C(:,2)'];
			depth =  [C(:,3)'];
			mag   =  [C(:,4)'];
			strike=  [C(:,5)'];
			dip   =  [C(:,6)'];
			rake  =  [C(:,7)'];
			size  =  [C(:,8)']; 
         end


         Focals.lola  = [long',lat'];
         Focals.depth = depth';
         Focals.mag   = mag';
		 Focals.strike= strike';
		 Focals.dip   = dip';
		 Focals.rake  =rake';
		 Focals.size  =size';
         Focals.xy    = [ lola2xy(Focals.lola',basemap,1) ]';
         ind=find(Focals.mag<0);    Focals.lola(ind,:)=[]; Focals.depth(ind)=[]; Focals.mag(ind)=[]; Focals.xy(ind,:)=[]; Focals.strike(ind,:)=[]; Focals.dip(ind,:)=[]; Focals.rake(ind,:)=[]; Focals.size(ind,:)=[]
         ind=find(Focals.depth==0); Focals.lola(ind,:)=[]; Focals.depth(ind)=[]; Focals.mag(ind)=[]; Focals.xy(ind,:)=[]; Focals.strike(ind,:)=[]; Focals.dip(ind,:)=[]; Focals.rake(ind,:)=[]; Focals.size(ind,:)=[]

     save(focals_matfile,'Focals')

