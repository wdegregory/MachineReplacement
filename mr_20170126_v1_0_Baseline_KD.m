lc
clear
dbstop if error

par.bigM=999999;
par.gamma=.95;
par.epsilon=.1;
par.numS=4;
par.numA=3;
par.alpha=.95;
N=2000;

%Transition Probability Matrices
data.TPM=[.7,.95,.999,1;
    0,.8,.98,1;
    0,0,.7,1;
    0,0,0,1];

%Cost Matrix
data.CM=[0,21000,40000;
    1300,21000,40000;
    3000,21000,40000;
    par.bigM,31000,40000];

Display=zeros(N,7);
Alpha=zeros(N,1);
S=zeros(N,1);
PDS=zeros(N,1);
V=zeros(N,par.numS);

tempV=zeros(1,par.numA);

S(2,1)=2;
Alpha(1,1)=par.alpha;

%Updated Value Vector
for n=2:N
    omega=rand(3);
    %Decision(1- do nothing, 2-return to state 1, 3-return to state 0)
    clear tempV;
    if S(n,1) == 4
        tempV(1,1:2)=par.bigM;
        tempPDS=1;
    else
        for a=1:par.numA-1
            %Make decision 2, return to state 1
            if a==1
                tempPDS=S(n,1);  
                %Make decision 3, return to state 0
            elseif a==2
                tempPDS=2;
                %No decision, stay in current state
            end
            tempV(1,a)=data.CM(S(n,1),a)+par.gamma*V(n-1,tempPDS);
        end
       tempPDS=1; 
    end
    tempV(1,par.numA)=data.CM(S(n,1),par.numA)+par.gamma*V(n-1,tempPDS);
    %Updated vector
    
    [vhat,aIndex]=min(tempV); 
    if aIndex == 3
        PDS(n,1) = 1;
    elseif aIndex == 2
        PDS(n,1) = 2;
    else
        PDS(n,1) = S(n,1);
    end
   
    %Updating through one step into the future, so a t is not necessary
    V(n,:)=V(n-1,:);
    if n>=3
        V(n,PDS(n-1,1))=(1-Alpha(n-1))*V(n-1,PDS(n-1,1))+Alpha(n-1)*vhat;
    end
    
    s=1;
    while s<=par.numS
        if omega(1)<=data.TPM(PDS(n,1),s)
            S(n+1)=s;
            s=par.numS+1;
        else
            s=s+1;
        end
    end
    
    Alpha(n)=max(Alpha(n-1)-1/N,0);
    
    Display(n,1)=S(n);
    Display(n,2)=aIndex; %action
    Display(n,3)=vhat;
    Display(n,4)=PDS(n);
    Display(n,5)=omega(1);
    Display(n,6)=S(n+1);
    Display(n,7)=Alpha(n);
    
    convCount=0;
    for s=1:par.numS
       if abs(V(n,s)-V(n-1,s))<par.epsilon 
           convCount=convCount+1;
       end
    end
    if n>10 && convCount==par.numS
        sprintf('converged at n = %d',n)
        break
    end
    
end
clear a r n s aIndex tempPDS tempV vhat N omega convCount