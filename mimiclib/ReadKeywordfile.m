function [S]=ReadKeywordfile(fname,delimiter);
%ReadKeywordfile     - reads an ascii keyword file 
%
%usage: [S] = ReadKeywordfile(fname,delimiter);
%
%Input:   fname      filename (*template, *rsc or *min)
%         delimiter  delimiter string (e.g. ' ', '=',':')   [default ' '] 
%
%Output:  S          structure containing keywords as fields
%
%    Read Keywords from an ascii file and puts contents into structure S. Reads *template and *min files (delimiter '=')
%    and *rsc files (delimiter ' '). Use the matlab isalpha function to determin wehter the value is a string or number. Exceptions 
%    are hardwired in the code (e.g. DATE12 and HEIGHT [23232-30000 is otherwise evaluated numerically). All strings containing '
%    are considered numbers so that they are evaluted as a proper matlab command (e.g.: strcat(getenv('GEODMODHOME'),'/DATADIR'))
%    Examples:
%
%   S=ReadKeywordfile('123.unw.rsc',' '):
%
%         ORBIT_NUMBER            38789-41190                   
%         VELOCITY                7552.29309370684              
%         HEIGHT                  0.7938263276E+06              
%         DATE                    030410
%
%   S=ReadKeywordfile('HawaiiRsatD1.template','='):
%
%         SLCdir                  = /RAID6/sbaker/SLCs/HawaiiRsatD6/
%         selectpairsopt.dir      = strcat(getenv('GEODMODHOME'),'/DATADIR')     
%         selectpairsopt.maxPbase = 600                    # some comment
%         selectpairsopt.minTbase = 1.0                             
%         electpairsopt.maxTbase  = 10    

%    Works with and # are comment lines.
%    Reads *template file 
%
%FA, Aug 2007    generalized from read_options_fromfile and Noel's ReadRscfileN.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~exist(fname,'file')    errordlg(['file ' fname ' does not exist']) ; error('user error -- exiting') ; end  
   if nargin~=2 delimiter=' '; end

%
% set keywords that are forced to be strings or numbers ignoring the isalpha command
%

   keyword_force_string_rsc      = {'DATE12' 'ORBIT_NUMBER'};                                 % to interprete e.g. 020112-031112 as string
   keyword_force_string_min      = {'makesaropt.igram2motionopt.tforstack' 'startDate'          'endDate'};
   keyword_force_string_min      = {keyword_force_string_min{:}            'masterPeriod'       'slavePeriod'};
   keyword_force_string_min      = {keyword_force_string_min{:}            'interferogramList'  'excludeList'};
   keyword_force_string_template = {''};

   keyword_force_number_rsc      = {'HEIGHT'};                                                % to interprete e.g. 0.7938263276E+06 as number
   keyword_force_number_min      = {''};
   keyword_force_number_template = {''};

   keyword_force_string          = [keyword_force_string_rsc keyword_force_string_min keyword_force_string_template]; 
   keyword_force_number          = [keyword_force_number_rsc keyword_force_number_min keyword_force_number_template]; 

   value_force_number            = {'true' 'false'}; 

%
% read file content into cell string
%
   fid=fopen(fname);   
   Call = textscan(fid,'%s','delimiter','\n'); 
   fclose(fid);

% convert into cell array containing the lines (don't know why this works) and replace #-comments by %-comments

   Clines = Call{1};                            % first convert into cell array containing the lines (% don't know why this works)
   Clines = strrep(Clines,'#','%');             

% remove comments and comment lines

   i=1;
   while i<=length(Clines)                     
         if any(strfind(Clines{i},'%')) Clines{i}=Clines{i}(1:strfind(Clines{i},'%')-1); end
         if isempty(Clines{i})  
            Clines(i)=[] ; 
         else
            i=i+1;
         end
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now put the remainder into optionsopt structure for return to calling function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for i=1:length(Clines)
      if length(Clines{i})                                                 % skip if empty line
         if isspace(delimiter)  
            C = textscan(strtrim(Clines{i}),'%s%s','commentStyle','%');    % this allows for multiple whitespace delimiter (not in second case)
         else 
            C = textscan(strtrim(Clines{i}),'%s%s','commentStyle','%','delimiter',delimiter);  
         end
         keyword = char(deblank(C{1}));
         value   = char(deblank(C{2}));
         if length(keyword) && length(value)                               % skip if keyword or value empty (don't do following if length(keyword)==0 || length(values)==0))
             %
             % find, evaluate and replace Unix environment variables (end of environment variable is blank or '/' (filesep), 
             % tested for value='$SLCDIR /RAID3/famelung', value='--geo --kmz --to $USER', value='$SLCDIR/Wellsquake'
             %
             if  any(strfind(value,'$'))
                   ind_dollar             = strfind(value,'$');

                   ind_filesep            = strfind(value(ind_dollar+1:end),filesep)-1;
                   ind_blank              = strfind(value(ind_dollar+1:end),' ')-1;
                   ind_end                = length(value(ind_dollar:end))-1;
                   environmentvar_length  = min([ind_blank ind_end ind_filesep]);
                   value_mod              = [value(1:ind_dollar-1) getenv(value(ind_dollar+1:ind_dollar+environmentvar_length)) value(ind_dollar+environmentvar_length+1:end) ];  

                   value                  = value_mod;
             end
             %
             % modify value field if keyword on force_string list
             %
             if     any(strcmp(keyword,keyword_force_string)) 
                        value_mod = sprintf('''%s''',value) ;              % force to string: put value into ' ' e.g. SLCdir = '/RAID6/sbaker/SLCs'
             elseif any(strcmp(keyword,keyword_force_number))  
                        value_mod = value;                                 % force to number: do not put into ' ', even if letter  (e.g. HEIGHT  0.7938263276E+06)
             elseif any(strcmp(value,value_force_number))  
                        value_mod = value;                                 % force to number: do not put into ' ', even if letter  (e.g. FollowGradient = true)
             elseif any(isalpha(value)) &&  ~any(strfind(value,'''')) 
                        value_mod = sprintf('''%s''',value) ;              % force all fields with letters to strings ' ', except if string contains ' (e.g. strcat(getenv('GEODMODHOME'),'/DIR')
             else
                         value_mod=value;
             end
             %disp(sprintf('S.%s=%s;',keyword,value_mod))  
             eval(sprintf('S.%s=%s;',keyword,value_mod));
         end
      end
   end


