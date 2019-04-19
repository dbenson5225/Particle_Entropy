clear all
L=12;  ttime=500;
nn=[30 100 300 1000 3000 10000 30000 ]
numsamp=15;
samples=L*(rand(numsamp,1)-0.5);
D=1e-3;
methodd='implicit'
Gau=@(x,mu,VAR) (1./sqrt(2.*pi.*VAR)).*exp(-((x-mu).^2)./(2.*VAR));

for k=1:length(nn)
    NN=nn(k)+1
    dx=L/(NN-1);  %dt=0.25*dx^2/D;  
    nt=0;
    A=zeros(NN,1);  active=ones(size(A)); active(1)=0; active(NN)=0;
    xtrue=(-L/2:dx:L/2)';    
    A(ceil(NN/2))=1/dx;
    time=0; dt=0.05;
    %  Figure out which FD cell is closest to the samples
     [idx r]=rangesearch(xtrue, samples, dx/2);
     for l=1:size(samples)         
         sampidx(l)=idx{l};
         dist(l)=r(l);
     end
    r=D*dt/dx^2; 
%    Amat = diag(-r*ones(NN-3,1),-1)+diag(1+2*r*ones(NN-2,1))+diag(-r*ones(NN-3,1),1);
%    Ainv=inv(Amat);     
        
%A=A.*active;
     
    while time<ttime;
        time=time+dt; nt=nt+1; tplot(nt)=time;
        [A] = diff_1d_spdiags(active,A,D,dt,dx,methodd);  
        %solve  
%        A(2)=A(2)+r*A(1);  A(NN-1)=A(NN-1)+r*A(NN);    
%        A(2:NN-1)=Ainv*A(2:NN-1);
       
        true=Gau(xtrue,0,2*D*time);
        true=true/(dx*sum(true));
 %       RMSEall(nt,k)=sqrt(mean((true-A).^2));
        SSEall(nt,k)=sum((true-A).^2);
        SSEsamp(nt,k)=sum((true(sampidx)-A(sampidx)).^2);
        corrSSE(nt,k)=sum(((true-A).^2)./true);
        dt=1.00*dt;
%    sum(A)*dx
%    sum(true)*dx
        
        figure(10)
        plot(xtrue,A,'o')
        hold on
        plot(xtrue,true,'-k');
        hold off
        drawnow
    end
    figure(11)
    plot(tplot,log(SSEsamp(:,k)/numsamp),'o');
    hold on 
    figure(21)
    plot(tplot,log(SSEall(:,k)/nn(k)),'sq');
    hold on

%    figure(11)
%    loglog(tplot,RMSEall(:,k),'o');
%    hold on
%    loglog(tplot,RMSEall(:,k)-log(dx),'o');
    
    figure(12)
    plot(tplot,-log(L/(nn(k)))+0.5*nn(k)*(1+log(2*pi)+nn(k)*log(SSEall(:,k)/nn(k))),'o');
    hold on
    figure(22)
    plot(tplot,0.5*nn(k)*(1+log(2*pi)+log(SSEall(:,k)/nn(k))),'sq');
    hold on
    figure(32)
    plot(tplot,0.5*nn(k)*(-log(L/nn(k))+1+log(2*pi)+log(SSEall(:,k)/nn(k))),'sq');
    hold on

    
    figure(13)
    plot(tplot,-log(L/(nn(k)))+numsamp*log(2*pi)+numsamp*log(SSEsamp(:,k)/numsamp)+numsamp*(SSEsamp(:,k)/numsamp),'o');
    hold on
    figure(23)
    plot(tplot,numsamp*log(2*pi)+numsamp*log(SSEsamp(:,k)/numsamp)+numsamp*(SSEsamp(:,k)/numsamp),'sq');
%    hold on

    %   loglog(tplot,SSEall(:,k),'o');
    
%    figure(13)
%    loglog(tplot,corrSSE(:,k));
%    hold on
%    loglog(tplot,corrSSE(:,k)-log(L/nn(k)),'o');

     drawnow
end
    %+log((1:100:nt)')
figure(1)
for k=1:length(nn)
    plot(tplot(1:100:end),-log(L/nn(k))+2*log(SSEall(1:100:end,k)/nn(k)),'o');
    hold on
end
legend('30','100','300','1000','3000','10000','30000')
figure(2)
for k=1:length(nn)
    plot(tplot(1:100:end),2*log(SSEall(1:100:end,k)/nn(k)),'o');
    hold on
end
legend('30','100','300','1000','3000','10000','30000')
    
figure(3)
for k=1:length(nn)
    %plot(tplot(1:100:end),-log(L/nn(k))+log((1:100:nt)')+log(SSEsamp(1:100:end,k)/numsamp),'o');
    plot(tplot(1:100:end),-log(L/nn(k))+2*log(SSEsamp(1:100:end,k)/numsamp),'o');
    hold on
end
figure(4)
for k=1:length(nn)
    plot(tplot(1:100:end),2*log(SSEsamp(1:100:end,k)/numsamp),'o');
    hold on
end

figure(5)
plot(log(nn),-log(L./nn)+2*log(SSEall(floor(end/2),:)./nn),'-o');
hold on 
plot(log(nn),2*log(SSEall(floor(end/2),:)./nn),'-sq');


