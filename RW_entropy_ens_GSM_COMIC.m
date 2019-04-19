clear all
%set(0, 'DefaultLineLineWidth', 2);
nn=[30 100 300 1000 3000 10000 30000 ];   % nn is a set of particle numbers to run
%nn=[300 3000];
mass=1; tottime=1000;  buffer=10;  everyother=1;   % buffer is "domain size"
ens=[80 80 60 50 20 6 4];     % Number of runs in each nn ensemble
%ens=[4 4];       

%sep=buffer/N;   binx=100*sqrt(dt)/N;  analdx=100*sqrt(dt)/N;
% H=-mp*N*log(mp);
D=1e-3;              % Diffusion coeff.
numinit=1;

eta=1;      %  fraction of D done by mass-transfer
entropy=zeros(300,5); comic=zeros(300,length(nn)+1);
SSEall=zeros(300,length(nn)); SSEvar=zeros(300,length(nn));

Gau=@(x,mu,VAR) (1./sqrt(2.*pi.*VAR)).*exp(-((x-mu).^2)./(2.*VAR));


for ll=1:length(nn)     % Different particle numbers loop
    ll
  for kens=1:ens(ll)    % Ensemble loop
      kens
    clear PRW; clear MT; clear Cmap;
    N=nn(ll);
    sep=buffer/N;   binx=buffer/N;  analdx=10*buffer/N;  mp=mass/N;
    PRW=zeros(N,1); MT=zeros(N,2);
    MT(1:N,1)=buffer*(rand(N,1)-0.5);
    MT(1:numinit,1)=0;  
    mp=1/numinit;
    MT(1:numinit,2)=mp;
    dt0=.01; dt=dt0;     % Initial delta t
    time=0; n=0; nt=0;
while time<tottime
   time=time+dt;   % I am allowing dt to grow over time (see bottom of time loop)
   n=n+1; nt=nt+1; tplot(nt)=time;   
%  Random Walk particles

if(eta<1)
    PRW(:,1)=PRW(:,1)+sqrt(2*D*dt)*randn(N,1);
%    MT(:,1) =MT(:,1)+sqrt(2*(1-eta)*D*dt)*randn(N,1);
    MT(:,1) =MT(:,1)+sqrt(2*(1-eta)*D*dt)*randn(N,1);
end
%  Mass transfer between particles
%  Deff is the effective diffusion per particle - this can included kernel
%  half-widths if desired.
    Deff=eta*D;     % Benson method (1-eta is RW, eta is MT)

    h_s=sqrt(4*Deff*dt);
    beta=1/2;
    
  % From a Sole-Mari paper; gsm: this is fine for interp.
    h_s = 0.822*((2*D*time)^(2/5))*(sep^(1/5));
    beta=2*Deff*dt/h_s^2;

    denom=2*h_s^2;
    factor= 1/sqrt(2*pi*h_s^2);
    
    dist=3*sqrt(denom/2);      % Max search radius for speed
%    factor=1/sqrt(8*pi*(Deff)*dt);denom=-8*(Deff)*dt;

    [idx r]=rangesearch(MT(:,1),MT(:,1),dist,'BucketSize',10);
  
  for i=1:N
     blah=idx{i};     % This is the indices list of nearby particles to i
     s=r{i};          % Associated radii
     Ptot=factor*exp(-s.^2/denom);
     rescale=sum(Ptot);  
     s=s(blah>i);
     Ptot=Ptot(blah>i);
     apple=blah(blah>i);
     if(length(apple)>0)
         order=randperm(length(apple));
         for j=1:length(apple)
            jpart=apple(order(j));    % Index of the current "other" particle
            v_s=Ptot(order(j))/rescale;
         
            dmA = beta*(MT(i,2)-MT(jpart,2))*v_s;
            %update particle masses of A 
            MT(i,2)=max(0,(MT(i,2)-dmA));
            MT(jpart,2)=max(0,(MT(jpart,2)+dmA));
        end   
     end  % If statement for more than one b particle 
  end  % Done with all particles

   entropy(n,1)=time; comic(n,1)=time; entropyMT(n,1)=time;
   maxx=5*sqrt(2*D*time);
%   binx=100*sqrt(time)/N;  analdx=100*sqrt(time)/N;
% Semi-analytic
   analdx=min(0.1,sqrt(2*D*time));
   xgau=-maxx:analdx:maxx; 
   truef=Gau(xgau,0,2*D*time);
   entropy(n,5)=-sum(analdx*truef.*log(truef));        % Inconsistent Entropy H_I
%   entropy(n,5)=-log(analdx)-sum(analdx*truef.*log(truef));        % Consistent Entropy H_C
%   entropy(n,5)=0.5*log(4*pi*D*time)+0.5;
%  Binned method
   ndx=ceil(2*maxx/binx);
   for j=1:ndx
       binx0=-maxx+(j-1)*binx;
       Cmap(j)=(1/N/binx)*sum(PRW(:,1)>binx0 & PRW(:,1)<binx0+binx);  
   end

%%%%%%%%%%  Entropies (OLD STYLE!!)
% Kernel method
   xkern=xgau;
   kernval=0*xkern;
% Use Gaussian with constant sigma
   varkern=2*D;
% Use Sole-Mari (2017)  
   varkern=(1.06*sqrt(2*D*time)*N^(-1/5))^2;


  for k=1:N
       kernval=kernval+Gau(xkern,PRW(k,1),varkern);
%       kernval=kernval+Gau(xkern,PRW(k,1),1*D*(n*dt)^.4);
  end
   kernval=(1/N)*kernval;                   % each particle has mass 1/N
   kernuse=kernval(kernval>1e-30);
   entropy4=-sum(analdx*kernuse.*log(kernuse));       % Inconsistent Entropy H_I
%   entropy4=-log(analdx)-sum(analdx*kernuse.*log(kernuse));       % Consistent Entropy H_C
   entropy(n,4)=((kens-1)*entropy(n,4)+entropy4)/kens;
%%%%%%%%%%
%  Random walk binned entropy 
   probPRW=binx*Cmap(Cmap>0);
   entropy2=-sum(probPRW.*log(probPRW));              % Consistent H_D
   entropy(n,2)=((kens-1)*entropy(n,2)+entropy2)/kens;  
%%%%%%%%%%
%  Mass-transfer alg. entropy
   probMT=MT(MT(:,2)>0,:);
   entropy3=-sum(probMT(:,2).*log(probMT(:,2)));      % Consistent H_D
   entropy(n,3)=((kens-1)*entropy(n,3)+entropy3)/kens;
   comic(n,ll+1)=entropy(n,3)+log(sep);
   entropyMT(n,ll+1)=entropy(n,3);
   
%%%%   Plot concentrations
 
%Not right now ...
%if(false)
%figure(1)
%xmapplot=-maxx:binx:-maxx+(ndx-1)*binx;
%   plot(xmapplot(1:everyother:end),Cmap(1:everyother:end),'gsq')
%   hold on 
%   plot(MT(1:everyother:end,1),MT(1:everyother:end,2)./sep,'rd')
%   plot(xkern(1:everyother:end),kernval(1:everyother:end),'bo-');
%   plot(xgau(1:everyother:end),truef(1:everyother:end),'-k','LineWidth',2);
%   hold off
%%%%%  Plot "dilution index".   
%   figure(2)
%   loglog(entropy(1:n,1),binx*exp(entropy(1:n,2)),'gsq')
%   hold on 
%   loglog(entropy(1:n,1),sep*exp(entropy (1:n,3)),'rd')
%   loglog(entropy(1:n,1),exp(entropy (1:n,4)),'bo')
%   loglog(entropy(1:n,1),sqrt(4*pi*D*exp(1)*entropy (1:n,1)),'-k')
%   %loglog(entropy(1:n,1),analdx*exp(entropy (1:n,5)),'-k')
%   axis([dt0 1000 .001 10])
%   hold off
%%%%%  Plot Entropies (WITH CONSISTENT CORRECTION).   
%   figure(3)
%   loglog(entropy(1:n,1),(entropy(1:n,2)),'gsq')   % Binned RW
%   hold on 
%   loglog(entropy(1:n,1),(entropy (1:n,3)),'rd')   % Mass-transfer
%   loglog(entropy(1:n,1),(entropy (1:n,4)),'bo')    % Kernel f(x)
%   loglog(entropy(1:n,1),(entropy (1:n,5)),'-k')    % Analytic f(x)

%   loglog(entropy(1:n,1),(log(N/buffer)+entropy (1:n,4)),'bo')    % Kernel f(x)
%   loglog(entropy(1:n,1),(log(N/buffer)+entropy (1:n,5)),'-k')    % Analytic f(x)
%   axis([dt0 1000 .1 10])
%   hold off

%   drawnow   
%end
   
%  Get SSEs for model fitness

   true=Gau(MT(:,1),0,2*D*time);
   true=true./sum(true)/sep;
   SSEnow=sum((true-MT(:,2)./sep).^2)/N;
   
   SSEall(nt,ll) = ((kens-1)*SSEall(nt,ll)+SSEnow)/kens;
   if kens>1
       SSEvar(nt,ll) =( (kens-2)*SSEvar(nt,ll) + (SSEall(nt,ll)-SSEnow).^2 ) /(kens-1);
   end 
   figure(90)
   plot(MT(:,1),MT(:,2)./sep,'rd');
   hold on
   plot(MT(:,1),true,'k+');
   legend('SSE/N = ',num2str(SSEall(nt,ll))); 
   hold off
   drawnow             
       
  dt=1.05*dt;
  
end  % timestep loop 

%  Plot dH/dt
figure(4)
loglog(0.5*(entropy(1:n-1)+entropy(2:n)),diff(entropy(1:n,2))./diff(entropy(1:n,1)),'gsq')
hold on
loglog(0.5*(entropy(1:n-1)+entropy(2:n)),diff(entropy(1:n,3))./diff(entropy(1:n,1)),'rd')
loglog(0.5*(entropy(1:n-1)+entropy(2:n)),diff(entropy(1:n,4))./diff(entropy(1:n,1)),'bo')
loglog(0.5*(entropy(1:n-1)+entropy(2:n)),diff(entropy(1:n,5))./diff(entropy(1:n,1)),'-k')
hold off

  end % end enemble loop
end   % end of different particle numbers
%plot ENTROPIES

figure(5)
for k=2:length(nn)+1
    plot(log(comic(1:1:n,1)),comic(1:1:n,k),'o');
    hold on
end

figure(6)
for k=2:length(nn)+1
    plot(log(comic(1:1:n,1)),entropyMT(1:1:n,k)+log(1/nn(k-1)),'sq');
    legend('100','300','1000','3000','10000','30000','Location','northwest')
    xlabel('ln(time)')
    ylabel('H-ln(n)')
    hold on
end
   plot(log(entropy(1:n,1)),(entropy (1:n,5))-log(buffer),'-k')
   hold off


figure(7)
for k=2:length(nn)+1
    plot(log10(comic(1:1:n,1)),comic(1:1:n,k),'sq');
    legend('300','1000','3000','10000','30000','Location','northwest')
    xlabel('log_{10}(time)')
    ylabel('H+ln(\Omega/n)')
    hold on
end
   plot(log10(entropy(1:n,1)),(entropy (1:n,5)),'-k')
   hold off

SSEall(nt+1:end,:)=[];   
SSEvar(nt+1:end,:)=[];
SSEstd=SSEvar.^0.5;
pos=log(SSEall+SSEstd)-log(SSEall);
neg=log(SSEall-SSEstd)-log(SSEall);

figure(21)
for k=1:length(nn)
    plot(tplot(1:end),-log(buffer/nn(k))+2*log(SSEall(1:end,k)),'o');
    hold on
end
legend('100','300','1000','3000','10,000','30,000','Location','northwest');
figure(22)
for k=1:length(nn)
    plot(tplot(1:end),2*log(SSEall(1:end,k)),'o');
    hold on
end
legend('100','300','1000','3000','10,000','30,000','Location','northwest');


figure(23)
errorbar(log(nn),-log(buffer./nn)+2*log(SSEall(end,:)),neg(end,:),pos(end,:),'-o')
hold on
errorbar(log(nn),2*log(SSEall(end,:)),neg(end,:),pos(end,:),'-o')



