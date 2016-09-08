function [fig_title,colorbar_title,ge_fname]=GenerateNamesFromDataStructure(data)
%GenerateNamesFromDataStructure   - generates name for plotdata plot, colorbar, googleearth
%
%usage: [fig_title,colorbar_title,ge_name]=GenerateNameFromDataStructure(data)
%
%Input:   data      name to which suffic will be appended  (e.g. rate_RsatA3 --> rate_RsatA3_qt)
%
%Output:  fig_name       used in PlotData
%         colorbar_title 
%         
%  Falk Amelung, June 2007
%
[DataSet,datestr,ratingstr,avgtimestr,stackstr,stackstr_ge,meanrating]=deal('');
Unit      = data.Unit;       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for title
           
   if isfield(data(1),'DataSet')             DataSet   = data(1).DataSet;                                        end
   if isfield(data(1),'date1')               datestr   = [data(1).date1(3:end) '-' data(1).date2(3:end) ];       end
   if isfield(data(1),'rating')              ratingstr = num2str(data(1).rating);                                end
   if isfield(data(1),'TotalTime')           avgtimestr= [' Averaged time:  ' num2str(data(1).TotalTime)      ]; end  
   if isfield(data(1),'StackTotalTime')      stackstr  = [' StackTotalTime: ' num2str(data(1).StackTotalTime) ];     
                                             stackstr_ge=['Stack'];                                              end
   if isfield(data(1),'MeanRating')          meanrating= [' MeanRating: ' num2str(data(1).MeanRating) ];         end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for colorbar
                                              prefix    ='LOS';
   if strcmp(DataSet,'Topography')            prefix    ='';                                                      end
   if strmatch(DataSet,{'East' 'North' 'Up'}) prefix    = 'ground';                                               end

                                              maintitle ='velocity';
   if strmatch(Unit,{'radian' 'm'})           maintitle ='displacement';                                          end
   if strcmp(DataSet,'Topography')            maintitle ='Topography';                                            end

                                              Unit      = data.Unit;       
   if isfield(data,'PlotUnit')                Unit      = data.PlotUnit;                                          end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig_title      = [ DataSet datestr ratingstr avgtimestr stackstr meanrating] ; 
    colorbar_title = [ prefix ' ' maintitle ' [' Unit ']' ] ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    list  = { DataSet datestr stackstr_ge} ; 

    fname ='';
    for i=1:length(list)
        if length(char(list(i)))  fname = [fname '_' char(list(i))] ; end
    end

    if length(fname)==0  
       ge_fname='GeodModImage';
    else
       ge_fname  = fname(2:end) ;    % cuts off the first '_'
    end

