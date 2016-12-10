function [projectname,areaname,satbeam,sat]=extract_ProjectName(name);
%extract_ProjectNamame,areaname,satbeam]names from string according to naming convention
%
%usage: [pname,aname,satbeam]=extract_ProjectName(fname);
%
%Input:   fname         directory name  (e.g. /RAID1/amelung/SO/HawaiiRsatA3_B0-1000-T300-400)
%
%Output:  ProjectName   project name  (e.g. HawaiiRsatA3)
%         areaname      name of the area (e.g. Hawaii)
%         satbeam       satellite viewing geometry (e.g. RsatA3)
%
%        examples: [projectname]=extract_ProjectName('/RAID6/amelung/HawaiiRsatD1') 
%                  [projectname]=extract_ProjectName('/RAID6/amelung/Stacks/POtest_HawaiiEnvD2_2002-2005') 
%
%        fname need to contain RsatF,RsatS,Rsat Env, Ers , Alos
%        e.g. RsatFA3, RsatD6, EnvA2, ErsA, Alos
%             (Radarsat standard beam is identified by Rsat or RsatS ) 
%        The program splits fname for '/', selects the last occurence containing Rsat,Env,Ers,Alos,Jers
%        and then splits the match for for '_'. If there is no match an empty string is returned.
%        This is used to extract some hardwired, beam related information (e.g. Vexcels Doppler Bandwidth in SBAS
%        and incidence angle in geodmod) 
%
%        Compare with: extract_ProjectName.pl
%
%  Falk Amelung, October 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pathstr, fname, ext] = fileparts(name);          % FA Dec 07: strip of extension
name = fullfile(pathstr,fname);

C = textscan(name,'%s','delimiter','/') ;         %split for '/'
    [Cf]=strfind(C{:},'Rsat'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Env') ; tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Ers' ); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Alos'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Jers'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Tsx') ; tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Csk') ; tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end %Anieri 5/16

if exist('match')
   name=match{end};
else
   projectname=[];areaname=[];satbeam=[];  % return emty string if there is no match
   return;
end

C = textscan(name,'%s','delimiter','_') ;        %split for '_'
    [Cf]=strfind(C{:},'Rsat'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Env') ; tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Ers' ); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Alos'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Jers'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Tsx'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end
    [Cf]=strfind(C{:},'Csk'); tmp=zeros(1,length(Cf)); for i=1:length(Cf) tmp(i)=~isempty(Cf{i}); end;  if find(tmp) match={C{1}{find(tmp)}};  end %Anieri 5/16


name=match{end};

str='Rsat' ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end
str='Env'  ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end
str='Ers'  ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end
str='Alos' ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end
str='Jers' ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end
str='Tsx'  ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end
str='Csk'  ; is=strfind(name,str) ;if (is) ie=is+length(str); projectname=name; areaname=name(1:is-1); satbeam=name(is:end); sat=str; end %Anieri 5/16
