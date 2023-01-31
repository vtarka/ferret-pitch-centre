% AUTHOR: Kerry Walker

function [MI,rMI,bMI,history]=MI2unbiased(h)

if isempty(h)
    MI=0;
else
    % matrix size h of the probability of each element
    ph=h/sum(h(:)); 

    % vector length = size(h,2), summed probability of each response (column)
    phx=sum(ph,1); 

    % vector length = size(h,1), summed probability of each stimulus (row)
    phy=sum(ph,2);

    % marginal distribution of stimuli and neural responses
    approx=phy*phx;

    % joint distribution of stimuli and neural responses
    sph=ph;
    sph(sph==0)=1;

    approx(approx==0)=1;

    % initial MI calculation
    rMI(1)=sum(sum(ph.*log2(sph./approx)));

    % initial bias calculation
    bMI(1)=bias(h);
    h2u=h;
    tmat{1}=h;
    step=1;
    while min(size(h2u))>1

        th=h2u;
        [mx,Ix]=min(phx); % minimum response probability
        [my,Iy]=min(phy); % minimum stimulus probability

        if my<mx % if the stimulus prob is less than the response prob

            th=th'; % transpose the matrix
            Ix=Iy; % flip the minimum indices (since we've transposed)
            phx=phy; % flip the stimuli and response probabilities

            % save the index of the minimum stimulus probability
            history(step,1)=2;
            history(step,2)=Iy;

        else
            history(step,1)=1;
            history(step,2)=Ix;
        end
        if Ix==1 % if the smallest stimulus probability is in the first column
            th(:,2)=th(:,2)+th(:,1); % add the first two columns together
            th=th(:,2:end); % eliminate the first column
        elseif Ix==length(phx) % if the smallest stim prob is the in the last column
            th(:,end-1)=th(:,end-1)+th(:,end); % add the last two columns together
            th=th(:,1:end-1); % eliminate the last column
        elseif phx(Ix-1)<phx(Ix+1)
            th(:,Ix-1)=th(:,Ix-1)+th(:,Ix);
            th=[th(:,1:Ix-1) th(:,Ix+1:end)];
        else
            th(:,Ix)=th(:,Ix)+th(:,Ix+1);
            th=[th(:,1:Ix) th(:,Ix+2:end)];
        end
        
        h2u=th;
        
        ph=h2u/sum(h2u(:)); % stim probability collapsed across responses
        phx=sum(ph,1); % response probability
        phy=sum(ph,2); % stim probability
        approx=phy*phx; % marginal probability
        sph=ph;
        sph(sph==0)=1;
        approx(approx==0)=1;
        rMI(end+1)=sum(sum(ph.*log2(sph./approx))); % new MI estimation
        bMI(end+1)=bias(h2u);
        tmat{end+1}=h2u;
    end
    tMI=rMI-bMI;
    MI=max(tMI(1:end-1));
end

function b=bias(h)
ss=size(h);
b=(ss(1)-1)*(ss(2)-1)/(2*sum(h(:))*log(2));
