  function [igram]=resize_igram(igram,opt)

%
% Modify the size of igram file
%
%
%  igram       :       igram structure
%
%   opt:
%
%       field		:		field to modify
%
%       subset      :       Subset igram structure 
%
%                       .ji or .lalo (indices and geographic coordinates respectively); 
%                     Can be: [y_top x_top y_length x_length]                             Top Left coordinate and size of box     
%                             [y_upleft x_upleft; y_lowright x_lowright]                  Top Left and Bottom Right coordinates
%                             [y_poly1 y_poly2 y_poly3 ...; x_poly1 x_poly2 x_poly3; ...] X and Y coordinates of polygone. Many polygones can be defined in which case
%                                              NaNs should placed in the list.
%
%                       Also accepts .i and .j and .la and .lo
%
%       resamp    	:       resamp igram structure opt.resamp=[xstep ystep]
%
%
% Example: opt = struct('field','data','subset',struct('lalo',[40 -120;38 -118]));
%          [out_igram] = resize_igram(igram,opt); 
%
%
% N. Gourmelen - June 2007
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(   'field'						   ,      {{'data'}} ,	   ...
                     'struct_test'                     ,        'off'    ,     ...
                     'subset' 	                       ,        'off'    ,     ...
                     'resamp'                      	   ,        'off'    )      ;


if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
for fi=1:length(field)
    
if ~isstruct(igram)
    struct_test=1;
    igram.data=igram;  igram.width=size(igram.data,1);  igram.file_length=size(igram.data,2);  
    igram.x_first=1;   igram.y_first=1;  igram.x_step=1; igram.y_step=1;
end

for ni=1:length(igram)

	Data = igram(ni).(genvarname(field{fi})); 
    
    if isstruct(subset)

        [XX,YY] = meshgrid(1:size(Data,2),1:size(Data,1)); 
        
        if isfield(subset,'ji')
            
            if size(subset.ji)==[1 4]
                
                [x1,y1] = meshgrid(subset.ji(2):subset.ji(2)+subset.ji(4)-1,subset.ji(1):subset.ji(1)+subset.ji(3)-1);
                points(:,1) = y1(:);  points(:,2) = x1(:);
                
            elseif size(subset.ji)==[2 2]
                
                [x1,y1] = meshgrid(subset.ji(2,1):subset.ji(2,2),subset.ji(1,1):subset.ji(1,2));
                points(:,1) = y1(:);  points(:,2) = x1(:);
                
            elseif size(subset.ji,1)==2
                
                in = inpolygon(XX(:),YY(:),subset.ji(2,:),subset.ji(1,:));  points = [YY(in),XX(in)];  
                
            else exit('Error, subset.ji has incorrect size')
                
            end

        elseif isfield(subset,'lalo')
            
            if size(subset.lalo)==[1 4]
                
                [upperleft] = LL2ind_igram(igram(ni),[subset.lalo(1) subset.lalo(2)]);
                [x1,y1]     = meshgrid(upperleft(2):upperleft(2)+subset.lalo(4)-1,upperleft(1):upperleft(1)+subset.lalo(3)-1);
                points(:,1) = y1(:);  points(:,2) = x1(:);
                
            elseif size(subset.lalo)==[2 2]
                
                [upleftlowdown] = LL2ind_igram(igram(ni),[subset.lalo(:,1) subset.lalo(:,2)]);
                [x1,y1]         = meshgrid(upleftlowdown(1,2):upleftlowdown(2,2),upleftlowdown(1,1):upleftlowdown(2,1));
                points(:,1)     = y1(:);  points(:,2) = x1(:);
                
            elseif size(subset.lalo,1)==2
                
                [indices] = LL2ind_igram(igram(ni),[subset.lalo(1,:)' subset.lalo(2,:)']);
                in = inpolygon(XX(:),YY(:),indices(:,2),indices(:,1));  points = [YY(in) XX(in)];
                
            else exit('Error, subset.ji has incorrect size')
                
            end

        elseif (isfield(subset,'j') && isfield(subset,'i'))
            
            points=[NaN NaN];
            for nj=1:size(subset.j,1)
                in=inpolygon(XX(:),YY(:),subset.i(nj,:)',subset.j(nj,:));  points=[points;YY(in) XX(in)];
            end
            points=points(2:size(points,1),:);
            
        elseif (isfield(subset,'la') && isfield(subset,'lo'))
            
            points=[NaN NaN];
            for nj=1:size(subset.la,1)
                [indices]=LL2ind_igram(igram(ni),[subset.la(nj,:)' subset.lo(nj,:)']);
                in=inpolygon(XX(:),YY(:),indices(:,2),indices(:,1));  points=[points;YY(in) XX(in)];
            end
            points=points(2:size(points,1),:);
            

        else error('Error, no proper subset field');

        end
        
        mask=repmat(NaN,size(Data,1),size(Data,2));  [L]=YX2L(points(:,1),points(:,2),size(mask,1));  mask(L)=1; 
        Data=Data(min(points(:,1)):max(points(:,1)),min(points(:,2)):max(points(:,2)),:);
        mask=mask(min(points(:,1)):max(points(:,1)),min(points(:,2)):max(points(:,2)));
        
        if length(size(Data))==3
            for nni=1:size(Data,3)
                
                Data(:,:,nni)=Data(:,:,nni).*mask;
            end
        else Data=Data.*mask;
        end

        igram(ni).(genvarname(field{fi}))=Data;  
        igram(ni).width=size(Data,2);igram(ni).file_length=size(Data,1);
        igram(ni).x_first=igram(ni).x_first+(min(points(:,2))-1)*igram(ni).x_step;
        igram(ni).y_first=igram(ni).y_first+(min(points(:,1))-1)*igram(ni).y_step;

    end

    if resamp
        resamp(2)=-abs(resamp(2));

        xfact=round(igram(ni).x_step/resamp(1)*size(igram(ni).(genvarname(field{fi})),2));yfact=round(igram(ni).y_step/resamp(2)*size(igram(ni).(genvarname(field{fi})),1));
        igram(ni).(genvarname(field{fi}))=resizem(igram(ni).(genvarname(field{fi})),[yfact xfact],'nearest');

        igram(ni).x_step=resamp(1);igram(ni).y_step=resamp(2);igram(ni).width=size(igram(ni).(genvarname(field{fi})),2);igram(ni).file_length=size(igram(ni).(genvarname(field{fi})),1);

    end

end

end % end length(field) loop

if struct_test  igram=igram.data;  end




