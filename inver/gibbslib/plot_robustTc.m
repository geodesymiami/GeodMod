function  plot_robustTc(models,energy,opt)
%calc_robustTc(FUN,bounds,optTc,varargin)
%
% FA may 2002, based on code from Peter Cervelli.

%Check OPTIONS
    if isempty(opt)
        error('calc_robust:  needs opt');
    else
        grid  =opt.igrid  ;
        matrix=opt.matrix;
        Tsched=opt.Tsched;
        runs  =opt.runs  ;
        %rstate=opt.rstate  ;
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
   strTc=sprintf('Tc: 10^%6.3f ', Tc );


% plot
number=size(mean_energy1,2)+1;
mpl=ceil(sqrt(number)); npl=ceil(number/mpl);

for r=1:number-1
  h=subplot(mpl,npl,r) ; 
  loglog(Tlist,mean_energy1(:,r)) ;  xlabel('Temperature'); ylabel('energy'); 
  title_str=sprintf('Run: %d',r);title(title_str)
%  semilogy(mean_energy1(:,r)) ; axis tight
  set(h,'XGrid','on','YGrid','on')
end
  h=subplot(mpl,npl,r+1) ;
  loglog(Tlist,mean_energy2) ;  xlabel('Temperature'); ylabel('mean energy')
  title_str=sprintf('Mean Energy over all runs: %d',r);title(title_str)
  title(['Mean Energy over all Runs  ' strTc])
  set(h,'XGrid','on','YGrid','on')
 
