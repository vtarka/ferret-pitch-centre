% AUTHOR: Kerry Walker

function [MI,rMI,bMI,history]=MI2unbiased(h)

if isempty(h)
    MI=0;
else
    ph=h/sum(h(:));
    phx=sum(ph,1);
    phy=sum(ph,2);
    approx=phy*phx;
    sph=ph;
    sph(sph==0)=1;
    approx(approx==0)=1;
    rMI(1)=sum(sum(ph.*log2(sph./approx)));
    bMI(1)=bias(h);
    h2u=h;
    tmat{1}=h;
    step=1;
    while min(size(h2u))>1
        th=h2u;
        [mx,Ix]=min(phx);
        [my,Iy]=min(phy);
        if my<mx
            th=th';
            Ix=Iy;
            phx=phy;
            history(step,1)=2;
            history(step,2)=Iy;
        else
            history(step,1)=1;
            history(step,2)=Ix;
        end
        if Ix==1
            th(:,2)=th(:,2)+th(:,1);
            th=th(:,2:end);
        elseif Ix==length(phx)
            th(:,end-1)=th(:,end-1)+th(:,end);
            th=th(:,1:end-1);
        elseif phx(Ix-1)<phx(Ix+1)
            th(:,Ix-1)=th(:,Ix-1)+th(:,Ix);
            th=[th(:,1:Ix-1) th(:,Ix+1:end)];
        else
            th(:,Ix)=th(:,Ix)+th(:,Ix+1);
            th=[th(:,1:Ix) th(:,Ix+2:end)];
        end
        
        h2u=th;
        
        ph=h2u/sum(h2u(:));
        phx=sum(ph,1);
        phy=sum(ph,2);
        approx=phy*phx;
        sph=ph;
        sph(sph==0)=1;
        approx(approx==0)=1;
        rMI(end+1)=sum(sum(ph.*log2(sph./approx)));
        bMI(end+1)=bias(h2u);
        tmat{end+1}=h2u;
    end
    tMI=rMI-bMI;
    MI=max(tMI(1:end-1));
end

function b=bias(h)
ss=size(h);
b=(ss(1)-1)*(ss(2)-1)/(2*sum(h(:))*log(2));
