function  [Tc,strTc]=find_robustTc(models,energy,opt)
%calc_robustTc(FUN,bounds,optTc,varargin)
%
% FA may 2002, based on code from Peter Cervelli.
% FA may 2007, quick+dirty update

%Check OPTIONS
    if isempty(opt)
        error('calc_robust:  needs opt');
    else
        grid  =opt.igrid  ;
        matrix=opt.matrix;
        Tsched=opt.Tsched;
        runs  =opt.runs  ;
        robustTc=opt.robustTc  ;
     end

   Tlist=fliplr(unique(Tsched));
   mean_energy1=zeros(length(Tlist),runs);
   mean_energy2=zeros(length(Tlist));
   for r=1:runs
       for i=1:length(Tlist)
           list=find(Tsched==Tlist(i));
           mean_energy1(i,r)=sum(energy(list,r)')/length(list) ;
       end
   end
   mean_energy2=mean(mean_energy1,2);
   i=find(mean_energy2(:)==min(mean_energy2(:)));
   Tc=log10(Tlist(i));

