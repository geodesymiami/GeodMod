function [Quakes]=MakeQuakes(quakefile,quakefileformat,basemap)
%MakeQuakes           - reads Quakefile and saves Quakes structure as matlab file
%
%  usage:  [Quakes]=MakeQuakes(quakefile,quakefileformat)
%
%  basemap is needed to switch from lola to xy
%
%TODO: should use quakeopt.file and quakeopt.format. Should introduce this the moment we work with 2 quakefiles.

global dir_out

  if quakefile(1) CheckInOut(quakefile,''); end             % check whether quakefile exist

  quakes_matfile = fullfile(dir_out,'quakes.mat');

  if  CheckInOut('',quakes_matfile) load(quakes_matfile); return; end
%
% open and read quakefile
%

  fid = fopen(quakefile);
  switch quakefileformat
  case{'anss_readable'}
        C=textscan(fid,'%s %s %f %f %f %f %s %s %s %s %s %s %s');
            date  = [C{1}'];
            time  = [C{2}'];
            lat   = [C{3}'];
            long  = [C{4}'];
            depth = [C{5}'];
            mag   = [C{6}'];            
                   % Notes:
                   %We need to read from the file according to the location (character) which is unchanged over the entire file.
                   %Also, how to specify header lines with textscan ?
                   %
                   %One way to to this would be:
                   %C=textscan(fid,'%s','delimiter','\n')
                   %tmp1= char(C{1}');       
                   %lat   = str2num(tmp1(:,25:30));
                   %
                   %but is this efficient ?
                   %
                   %Some information may be given my matlab's importwizard.
         case{'anss'}
            C=textscan(fid,'%f%f%f%s%f%f%s%s%s%s%s%s%s%s%s')
            year  = [C{1}'];
            month = [C{2}'];
            day   = [C{3}'];
            tmp1= char(C{4}');
            tmp2= char(C{5}');         % This does not work. Don't know why. Not necessary to fix if we use anss_readable
            mag   = [C{6}'];
            lat   = str2num(tmp1(:,9:15));
            %long  = str2num(tmp2(:,1:7));
            depth = str2num(tmp2(:,9:13));
         case{'INGV-CT-FocMecs'}
            C=textscan(fid,'%d %f %s %f %f %f %s %s %s %s %s %s','headerLines',2);
            date  = [C{2}'];
            time  = [C{3}'];
            lat   = [C{4}'];
            long  = [C{5}'];
            depth = [C{6}'];
            mag   = depth*0+1 ;    % The file I got from Amalia did not have an entry for magnitude
         case{'UNR'}
            C=textscan(fid,'%f %f %f %f %d %d %s %s','headerLines',0);
            lat   = [C{1}'];
            long  = [C{2}'];
            depth = [C{3}'];
            mag   = [C{4}'];
            date  = [C{7}'];
            time  = [C{8}'];
            %for i=1:length(date) datestring{i}=[date{i} ' ' time{i}]; end;
            %qq=datestr(datestring,29);
         case{'UNR2'}
            fclose all;
            C     =  xlsread(quakefile,1,'range','basic');
            lat   =  C(:,2)';
            long  =  C(:,1)';
            depth = -C(:,3)/1000';
            mag   =  C(:,4)';
         case{'RSMAS-SHORT'}
            C     =  textscan(fid,'%f %f %f %f','headerLines',0);
            long  =  [C{1}'];
            lat   =  [C{2}'];
            depth =  [C{3}'];
            mag   =  [C{4}'];
         end


         Quakes.lola  = [long',lat'];
         Quakes.depth = depth';
         Quakes.mag   = mag';
         Quakes.xy    = [ lola2xy(Quakes.lola',basemap,1) ]';
         ind=find(Quakes.mag<0);    Quakes.lola(ind,:)=[]; Quakes.depth(ind)=[]; Quakes.mag(ind)=[]; Quakes.xy(ind,:)=[];
         ind=find(Quakes.depth==0); Quakes.lola(ind,:)=[]; Quakes.depth(ind)=[]; Quakes.mag(ind)=[]; Quakes.xy(ind,:)=[];

         % 10/2008: the following is a perferred format but gave errors in
         % process_defaultoption
         %Quakes = struct('lola',   num2cell([long ;lat]',2)',   ...
         %                'mag',  num2cell(mag),               ...
         %               'depth',num2cell(depth));           
         %  ind = find([Quakes(:).mag]  <0); Quakes(ind)=[];
         %  ind = find([Quakes(:).depth]<0); Quakes(ind)=[]; 

     save(quakes_matfile,'Quakes')

