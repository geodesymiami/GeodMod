function [amp,phase,RscVar] = readfile(infile,opt)

%
%   readfile        - reads roi_pac, doris, gamma, IREA, ... file types
%
% usage:  [a,p,rscinfo] = readfile(infile,opt);
%
% opt can be:
%
%    - subset:
%         extracts subset if given:
%                  subset=[ymin xmin ypix xpix]
%         for subset=struct('YXLW',[100 100 200 300],'Proj','CART or GEO')  data(100:300,100:400) are extracted
%
%    - origin:
%	  Can be roi_pac (default), gamma, doris, ...
%
%    - dimension
%
% NG, August 07
%

defaultopt=struct(                                    ...
    'dimension'      ,        'off'      ,            ...
    'subset'         ,        'off'      ,            ...
    'precision'      ,        'double'   ,            ...
    'origin'         ,        'off'      )             ;

if exist('opt')
    [opt]=process_defaultoptions(opt,defaultopt);  display(opt)
else
    [opt]=process_defaultoptions('',defaultopt);  display(opt)
end
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end

%%%%%%%%%%%%%%%%%%%%

if ~origin
%    if (length(strfind(infile,'x')) && length(strfind(infile,'dat')) && length(strfind(infile,'_')))
%        origin='irea';
%    else
        origin='roi_pac';
%    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(origin)

    case{'gamma'}
        
        if ~dimension
            exit('Error, you need to provide dimensions');
        end
        infile=deblank(infile);fid=fopen(infile,'r','ieee-be');
        [F,g]=fread(fid,'float','ieee-be');
        indx=[1:2:size(F,1)];indy=[2:2:size(F,1)];
        amp=sqrt(F(indx).*F(indx)+F(indx).*F(indx));
        pha=atan2(F(indy),F(indx));
        pha=reshape(pha,dimension(2),dimension(1));
        amp=reshape(amp,dimension(2),dimension(1));
        phase=flipud(rot90(pha));amp=flipud(rot90(amp));
        
        
    case{'irea'}
        
        [RscVar]=Read_IREA(infile);  phase=RscVar.data;
        amp=repmat(NaN,size(phase));  RscVar=rmfield(RscVar,'data');

    case{'roi_pac'}

        infile=deblank(infile);extens=infile(size(infile,2)-2:size(infile,2));
        RscVar=ReadKeywordfile([infile,'.rsc']);

        if (strcmp(extens,'cor') || strcmp(extens,'unw') || strcmp(extens,'hgt'))

            switch(precision)
                case('double')
                    fid=fopen(infile,'r');
                    [F,count]=fread(fid,'float32');
                    fclose(fid);
                    
                    nlines   = length(F)/(2*RscVar.WIDTH);
                    recl     = 8 * RscVar.WIDTH;
                    
                    phase = zeros(nlines,RscVar.WIDTH);
                    amp = zeros(nlines,RscVar.WIDTH);
                    
                    for line = 1:nlines
                        amp(line,:)  = F((line-1)*(2*RscVar.WIDTH)+1:(line-1)*(2*RscVar.WIDTH)+RscVar.WIDTH)';
                        phase(line,:)= F((line-1)*(2*RscVar.WIDTH)+RscVar.WIDTH+1 : line*(RscVar.WIDTH*2))';
                    end
                case('single')
                    fid=fopen(infile,'r');
                    [F,count]=fread(fid,'float32=>single');
                    fclose(fid);
                    
                    nlines   = length(F)/(2*RscVar.WIDTH);
                    recl     = 8 * RscVar.WIDTH;
                    
                    phase = single(zeros(nlines,RscVar.WIDTH));
                    amp = single(zeros(nlines,RscVar.WIDTH));
                    
                    for line = 1:nlines
                        amp(line,:)  = F((line-1)*(2*RscVar.WIDTH)+1:(line-1)*(2*RscVar.WIDTH)+RscVar.WIDTH)';
                        phase(line,:)= F((line-1)*(2*RscVar.WIDTH)+RscVar.WIDTH+1 : line*(RscVar.WIDTH*2))';
                    end
            end
            

        elseif (strcmp(extens,'int') || strcmp(extens,'slc'))

            fid=fopen(infile,'r');
            [F,count]=fread(fid,'float');

            indx=[1:2:size(F,1)];indy=[2:2:size(F,1)];
            amp=sqrt(F(indx).*F(indx)+F(indy).*F(indy));
            pha=atan2(F(indy),F(indx));
            pha=reshape(pha,RscVar.WIDTH,RscVar.FILE_LENGTH);
            amp=reshape(amp,RscVar.WIDTH,RscVar.FILE_LENGTH);
            phase=flipud(rot90(pha));amp=flipud(rot90(amp));


        elseif (strcmp(extens,'dem'))

            fid=fopen(infile,'r');
            [dem]=fread(fid,[RscVar.WIDTH RscVar.FILE_LENGTH],'short');
            demflr=fliplr(dem);phase=rot90(demflr);amp=phase*0;

        else

            error(sprintf('%s: exiting --- Unknown file type'));

        end

end

if isstruct(subset)

    [x1,y1]=meshgrid(1:size(phase,2),1:size(phase,1)) ;

    subset.x=[subset.x subset.x(:,1)];subset.y=[subset.y subset.y(:,1)];  points=[NaN,NaN];

    if isfield(subset,'Proj')==0

        for ni=1:size(subset.x,1)
            in=inpolygon(x1(:),y1(:),subset.x(ni,:),subset.y(ni,:));  points=[points(:,1),points(:,2);y1(in),x1(in)];
        end

        RscVar.X_FIRST=RscVar.X_FIRST+(min(points(:,2))-1)*RscVar.X_STEP;RscVar.Y_FIRST=RscVar.Y_FIRST+(min(points(:,2))-1)*RscVar.Y_STEP;

    elseif strcmp(subset.Proj,'GEO')==1

        for ni=1:size(subset.x,1)

            locind=LL2ind_hgt(infile,[subset.y(ni,:)' subset.x(ni,:)']);
            in=inpolygon(x1(:),y1(:),locind(:,2),locind(:,1));  points=[points(:,1),points(:,2);y1(in),x1(in)];
            xall(:,ni)=locind(:,2);yall(:,ni)=locind(:,1);

        end

        junk=subset.x(find(xall'==min(min(xall))));  RscVar.X_FIRST=junk(1);
        junk=subset.y(find(yall'==min(min(yall))));  RscVar.Y_FIRST=junk(1);

    else

        for ni=1:size(subset.x,1)
            in=inpolygon(x1(:),y1(:),subset.x(ni,:),subset.y(ni,:));  points=[points(:,1),points(:,2);y1(in),x1(in)];
        end

        RscVar.X_FIRST=RscVar.X_FIRST+(min(points(:,2))-1)*RscVar.X_STEP;RscVar.Y_FIRST=RscVar.Y_FIRST+(min(points(:,2))-1)*RscVar.Y_STEP;

    end

    points(1,:)=[];

    L=YX2L(points(:,1),points(:,2),size(phase,1));hh=1:size(phase,1)*size(phase,2);hh(L)=[];

    phase(hh)=NaN; phase=phase(min(points(:,1)):max(points(:,1)),min(points(:,2)):max(points(:,2)));
    amp(hh)=NaN; amp=amp(min(points(:,1)):max(points(:,1)),min(points(:,2)):max(points(:,2)));

    RscVar.WIDTH=size(phase,2);RscVar.FILE_LENGTH=size(phase,1);

end





