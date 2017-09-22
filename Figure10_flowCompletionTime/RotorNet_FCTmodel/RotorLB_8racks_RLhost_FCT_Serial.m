function []=RotorLB_8racks_RLhost_FCT_Serial()

% ------------------------------------------------------------------------%
% Copyright 2017 Regents of the University of California

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:

% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.

% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived
%    from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ------------------------------------------------------------------------%

%% Inputs:

fs=.2e6; % flow size, bytes

Nrack=8; % total # racks
Nhost=16; % total # hosts
Nhpr=Nhost/Nrack; % # hosts / rack

L=10e9/8; % link rate, Bytes / seconds
slot=200e-6; % slot time, seconds
duty=.9; % duty cycle: up time / total time

linkcap=int64(duty*slot*L); % bytes / ToR uplink in 1 slot
hostcap=int64(duty*slot*L); % bytes / host uplink in 1 slot

Ncycles=100;

Nflows=2240;

%% set up matchings

% assign number of switches (NS) and number of matchings per switch (Nmatch):

Nmatch=4; % number of matchings per switch (crossover matchings)
NS=ceil(Nrack/Nmatch); % number of switches

% construct bidirectional matchings
mc=cell(1,log2(Nrack));
for a=1:log2(Nrack)
    mc{a}=zeros(Nrack);
    for b=1:Nrack/(2^a)
        for c=0:(2^a-1)
            mc{a}((2^a)*b-c,(2^a)*b-((2^a-1)-c))=1;
        end;
    end;
end;

mc{4}=mc{3};
mc{3}=zeros(Nrack);
mc{3}(1,6)=1;
mc{3}(2,5)=1;
mc{3}(3,8)=1;
mc{3}(4,7)=1;
mc{3}(5,2)=1;
mc{3}(6,1)=1;
mc{3}(7,4)=1;
mc{3}(8,3)=1;

m=cell(1,Nmatch); % all connections seen at each of Nmatch timeslots
for a=1:Nmatch
    m{a}=mc{a};
    m{a}=m{a}+fliplr(eye(Nrack))*mc{a};
    % ensure main diagonal is zero
    m{a}(logical(diag(ones(1,Nrack),0)))=zeros(Nrack,1);
end;

%% bijective host-host demand:

% make sure hosts don't send to themselves or within the same rack
diag0=~diag(ones(1,Nhost));
for a=1:Nrack
    hostids=(a-1)*Nhpr+1:(a-1)*Nhpr+Nhpr;
    diag0(hostids(1),hostids(2))=0;
    diag0(hostids(2),hostids(1))=0;
end;

Dhost=diag0;

[src0,dst0]=find(Dhost);

%% flows

flows.src=repmat(src0.',1,10);
flows.dest=repmat(dst0.',1,10);
flows.strt=199.999e-6*ones(1,Nflows); % time, seconds 
flows.size=int64(fs*ones(1,Nflows)); % bytes

names={'src','dest','strt','size'};

% order flows by start time
[~,inds]=sort(flows.strt);
for a=1:4
    flows.(names{a})=flows.(names{a})(inds);
end;

% remove any flows that start after the simulation ends
if max(size(flows.src))>1
    indcut=find(flows.strt>(Ncycles*Nmatch)*slot,1);
    if isempty(indcut)==0
        for a=1:4
            flows.(names{a})=flows.(names{a})(1:indcut-1);
        end;
    end;
end;

% find src-dest pair with the most flows
maxflows=0;
usrc=unique(flows.src);
[~,nus]=size(usrc); % number of unique sources
for a=1:nus
    dests=flows.dest(flows.src==usrc(a));
    udest=unique(dests);
    [~,nud]=size(udest);
    for b=1:nud
        temp=sum(dests==udest(b));
        maxflows=max([maxflows temp]);
    end;
end;

% allocate 3D arrays to hold flow size and arrival times
flow_size=int64(zeros(Nhost,Nhost,maxflows));
flow_arrival=nan(Nhost,Nhost,maxflows);

for a=1:nus
    vindsrc=flows.src==usrc(a);
    dests=flows.dest(vindsrc);
    udest=unique(dests);
    [~,nud]=size(udest);
    srcinds=find(vindsrc);
    for b=1:nud
        vinddest=dests==udest(b);
        nelem=sum(vinddest);
        flow_size(usrc(a),udest(b),1:nelem)=flows.size(srcinds(vinddest));
        flow_arrival(usrc(a),udest(b),1:nelem)=flows.strt(srcinds(vinddest));
    end;
end;

% generate input demand per slot
slotcnt=0;
D=cell(Ncycles,Nmatch);
for a=1:Ncycles
    for b=1:Nmatch
        slotcnt=slotcnt+1;
        D{a,b}=flow_size.*int64((flow_arrival>=(slotcnt-1)*slot & flow_arrival<slotcnt*slot));
    end;
end;


%%

% initialize queues - index each flow by (src host, dest host, flow #)

Dworking=int64(zeros(Nhost,Nhost,maxflows)); % working demand
Rworking=int64(zeros(Nhost,Nhost,maxflows)); % working received data
flow_completion=nan(Nhost,Nhost,maxflows); % seconds

queues=cell(1,Nhost);
for a=1:Nhost
    queues{a}=int64(zeros(Nhost,Nhost,maxflows)); % indirect queues
end;

%% Simulator:

mscurrent=0; % current matching slot

delivered=zeros(Nmatch,Ncycles); % aggregate # bytes delivered / slot

for a=1:Ncycles % iterate through cycles of matchings
    
    fprintf('\nCycle %d',a);
    
    for b=1:Nmatch % iterate through matching slots
        
%         fprintf('\nCycle %d, slot %d',a,b);
        
        fine_latency=zeros(Nhost,Nhost,maxflows);
        
        Dworking=Dworking+D{a,b}; % current accrued demand
        
        mscurrent=mscurrent+1; % update matching slot
        
        % initialize link capacities for this slot
        con=find(m{b}(1,:)>0);
        [~,Ncon]=size(con); % get number of connections
        R=linkcap*int64(ones(Nrack,Ncon)); % bytes / slot, ToR uplinks
        Rhost=hostcap*int64(ones(Nrack,Nhpr)); % bytes / slot, host uplinks
        
        % ---------- phase 1 ---------- %
        
        % first, serve any existing nonlocal (two hop) traffic
        % abstract away where this data is being pulled from
        % ! assume this doesn't interfere with hosts sending traffic later
        for c=1:Nrack % sweep racks
            con=find(m{b}(c,:)>0); % get ids of connected racks
            for d=1:Ncon % sweep connected racks
                
                dhostids=(con(d)-1)*Nhpr+1:(con(d)-1)*Nhpr+Nhpr; % hosts in connected rack
                
                [queuesent,queues{c}(:,dhostids,:),R(c,d)]=fifo(queues{c}(:,dhostids,:),R(c,d));
                
                if sum(sum(sum(queues{c}(:,dhostids,:),3),2))~=0
                    fprintf('\nError! Nonlocal queues did not drain\n');
                    fprintf('sent = %d, remain = %d\n',sum(sum(sum(queuesent,3),2)), ...
                        sum(sum(sum(queues{c}(:,dhostids,:),3),2)));
                    fprintf('a = %d, b = %d, c = %d, con = %d',a,b,c,con(d));
                    return;
                end;
                
                if sum(sum(sum(queuesent,3),2))>0 % something was sent
                    % check if the flow(s) completed
                    Rworking(:,dhostids,:)=Rworking(:,dhostids,:)+queuesent;
                    inds=find(queuesent>0); % get linear indices
                    [i,j,k]=ind2sub([Nhost,Nhpr,maxflows],inds); % get matrix subscripts
                    [nsubs,~]=size(i);
                    addlatency=0; % additional latency due to serial delivery
                    for e=1:nsubs
                        addlatency=addlatency+double(queuesent(i(e),j(e),k(e)))/double(linkcap);
                        if addlatency>fine_latency(i(e),(con(d)-1)*Nhpr+j(e),k(e))
                            fine_latency(i(e),(con(d)-1)*Nhpr+j(e),k(e))=addlatency;
                        end;
                        if Rworking(i(e),(con(d)-1)*Nhpr+j(e),k(e))==flow_size(i(e),(con(d)-1)*Nhpr+j(e),k(e))
                            flow_completion(i(e),(con(d)-1)*Nhpr+j(e),k(e))= ...
                                slot*(mscurrent+fine_latency(i(e),(con(d)-1)*Nhpr+j(e),k(e)));
                        end;
                    end;
                    
                    delivered(b,a)=delivered(b,a)+sum(sum(sum(queuesent,3),2));
                    
                end;
            end;
        end;
        
        
        % next, send local data directly if possible
        % !!! account for hosts' sending capacities!
        for c=1:Nrack % sweep racks
            shostids=(c-1)*Nhpr+1:(c-1)*Nhpr+Nhpr; % sending host ids
            con=find(m{b}(c,:)>0); % get ids of connected racks
            
            % step 1: initial fairshare
            Rhostalloc=int64(zeros(Nhpr,Ncon));
            for d=1:Ncon
                dhostids=(con(d)-1)*Nhpr+1:(con(d)-1)*Nhpr+Nhpr;
                for e=1:Nhpr
                    temp=sum(sum(Dworking(shostids(e),dhostids,:),3),2);
                    Rhostalloc(e,d)=linkcap*int64(temp>0);
                end;
                % fairshare ToR uplink capacity across hosts
                Rhostalloc(:,d)=fairshare1(Rhostalloc(:,d),R(c,d));
            end;
            for d=1:Nhpr
                % fairshare host uplink capacity across ToR uplinks
                Rhostalloc(d,:)=fairshare1(Rhostalloc(d,:),Rhost(c,d));
            end;
            
            chk=1;
            while chk==1
                
                Rhostalloc0=Rhostalloc; % record the allocation before sending
                
                % step 2: send traffic according to Rhostalloc
                for d=1:Ncon
                    if R(c,d)>0 % there is remaining capacity on the link
                        % additional latency into the slot (may have already sent traffic)
                        addlatency=double(linkcap-R(c,d))/double(linkcap);
                        
                        dhostids=(con(d)-1)*Nhpr+1:(con(d)-1)*Nhpr+Nhpr;
                        
                        for e=1:Nhpr
                            [Dsent,Dworking(shostids(e),dhostids,:),Rhostalloc(e,d)]= ...
                                fifo(Dworking(shostids(e),dhostids,:),Rhostalloc(e,d));
                            
                            if sum(sum(Dsent,3),2)>0 % something was sent
                                % check if the flow(s) completed
                                Rworking(shostids(e),dhostids,:)=Rworking(shostids(e),dhostids,:)+Dsent;
                                inds=find(Dsent>0); % get linear indices
                                [~,j,k]=ind2sub([1,Nhpr,maxflows],inds); % get matrix subscripts
                                [nsubs,~]=size(j);
                                for f=1:nsubs
                                    addlatency=addlatency+double(Dsent(1,j(f),k(f)))/double(linkcap);
                                    if addlatency>fine_latency(shostids(e),(con(d)-1)*Nhpr+j(f),k(f))
                                        fine_latency(shostids(e),(con(d)-1)*Nhpr+j(f),k(f))=addlatency;
                                    end;
                                    if Rworking(shostids(e),(con(d)-1)*Nhpr+j(f),k(f))==flow_size(shostids(e),(con(d)-1)*Nhpr+j(f),k(f))
                                        flow_completion(shostids(e),(con(d)-1)*Nhpr+j(f),k(f))= ...
                                            slot*(mscurrent+fine_latency(shostids(e),(con(d)-1)*Nhpr+j(f),k(f)));
                                    end;
                                end;
                                
                                delivered(b,a)=delivered(b,a)+sum(Dsent(:));
                                
                            end;
                        end;
                    end;
                end;
                
                sent=Rhostalloc0-Rhostalloc;
                
                % step 3: update fairshare and check for completion
                
                for d=1:Ncon % update R
                    R(c,d)=R(c,d)-int64(sum(sent(:,d)));
                end;
                for d=1:Nhpr % update Rhost
                    Rhost(c,d)=Rhost(c,d)-int64(sum(sent(d,:),2));
                end;
                
                % prevent runaway loop:
                % fairshare1 is conservative, so there may be some
                % remaining capacity that can't be allocated to demand
                mask=int64(sent>0); % elements that used allocated capacity
                
                Rhostalloc=int64(zeros(Nhpr,Ncon));
                for d=1:Ncon
                    dhostids=(con(d)-1)*Nhpr+1:(con(d)-1)*Nhpr+Nhpr;
                    for e=1:Nhpr
                        temp=sum(sum(Dworking(shostids(e),dhostids,:),3),2);
                        Rhostalloc(e,d)=linkcap*int64(temp>0);
                    end;
                end;
                
                % apply mask
                Rhostalloc=Rhostalloc.*mask;
                
                for d=1:Ncon
                    % fairshare ToR uplink capacity across hosts
                    Rhostalloc(:,d)=fairshare1(Rhostalloc(:,d),R(c,d));
                end;
                for d=1:Nhpr
                    % fairshare host uplink capacity across ToR uplinks
                    Rhostalloc(d,:)=fairshare1(Rhostalloc(d,:),Rhost(c,d));
                end;
                
                % check whether to exit loop
                if sum(Rhostalloc(:))==0
                    chk=0;
                end;
                
            end;
            
        end;
        
        % ---------- phase 2 ---------- %
        % !!! account for hosts' sending capacities!
        
        allalloc=cell(1,Nrack); % store allocation arrays by sending rack
        % each cell has Ncon entries
        allcon=cell(1,Nrack); % store which rack is allocating
        % each cell has Ncon entries
        
        % bidirectional step - each rack computes how much indirect
        % data it can accept from each sending rack, given its own queues and the senders' offers
        for c=1:Nrack % sweep racks
            
            con=find(m{b}(c,:)>0); % get ids of connected racks
            % get local queues of connected nodes
            offer=int64(zeros(Ncon*Nhpr,Nhost));
            for d=1:Ncon
                offer((d-1)*Nhpr+1:(d-1)*Nhpr+Nhpr,:)= ...
                    int64(sum(Dworking((con(d)-1)*Nhpr+1:(con(d)-1)*Nhpr+Nhpr,:,:),3));
            end;
            
            % calculate available space at local rack
            bufferavail=int64(zeros(1,Nrack));
            for d=1:Nrack
                % account of existing local traffic:
                for e=1:Nhpr
                    temp=int64(sum(sum(Dworking((c-1)*Nhpr+e,(d-1)*Nhpr+1:(d-1)*Nhpr+Nhpr,:),2),3));
                    bufferavail(d)=bufferavail(d)+min([temp hostcap]);
                end;
                % account for existing nonlocal traffic:
                bufferavail(d)=bufferavail(d)+ ...
                    int64(sum(sum(sum(queues{c}(:,(d-1)*Nhpr+1:(d-1)*Nhpr+Nhpr,:),3),2)));
            end;
            bufferavail=linkcap*int64(ones(1,Nrack))-bufferavail;
            bufferavail(bufferavail<0)=0;
            bufferavail(c)=0;
            
            inds=find(bufferavail==0); % remove any offers where bufferavail = 0
            [~,Ninds]=size(inds);
            for d=1:Ninds
                offer(:,(inds(d)-1)*Nhpr+1:(inds(d)-1)*Nhpr+Nhpr)= ...
                    int64(zeros(Ncon*Nhpr,Nhpr));
            end;
            
            offerscl=int64(zeros(Ncon*Nhpr,Nhost));
            
            Rhostfs=int64(zeros(Ncon,Nhpr));
            for d=1:Ncon
                Rhostfs(d,:)=Rhost(con(d),:);
            end;
            for d=1:Ncon
                indtemp=find(m{b}(con(d),:)>0);
                Rhostfs(d,:)=fairshare1(Rhostfs(d,:),R(con(d),indtemp==c));
                for e=1:Nhpr
                    offerscl((d-1)*Nhpr+e,:)= ...
                        fairshare1(offer((d-1)*Nhpr+e,:),Rhostfs(d,e));
                end;
            end;
            
            alloc=int64(zeros(Ncon*Nhpr,Nhost)); % bandwidth allocated to senders
            chk=1;
            cnt=0;
            while chk==1
                cnt=cnt+1;
                if cnt==100
                    return;
                end;
                
                aggtemp=int64(zeros(1,Nrack));
                for d=1:Nrack
                    aggtemp(d)=int64(sum(sum(offerscl(:,(d-1)*Nhpr+1:(d-1)*Nhpr+Nhpr),2)));
                end;
                slices=find(aggtemp>0); % which groups to consider fairsharing over
                [~,Nslices]=size(slices);
                if Nslices>0
                    for d=1:Nslices
                        intemp=offerscl(:,(slices(d)-1)*Nhpr+1:(slices(d)-1)*Nhpr+Nhpr);
                        outtemp=fairshare1(intemp(:),bufferavail(slices(d)));
                        bufferavail(slices(d))=bufferavail(slices(d))-sum(outtemp);
                        if sum(outtemp)==0
                            bufferavail(slices(d))=0; % prevent runaway loop
                        end;
                        alloc(:,(slices(d)-1)*Nhpr+1:(slices(d)-1)*Nhpr+Nhpr)= ...
                            alloc(:,(slices(d)-1)*Nhpr+1:(slices(d)-1)*Nhpr+Nhpr)+ ...
                            reshape(outtemp,[Ncon*Nhpr Nhpr]);
                    end;
                    offerscl=offer-alloc;
                    
                    inds=find(bufferavail==0); % remove any offers where bufferavail = 0
                    [~,Ninds]=size(inds);
                    for d=1:Ninds
                        offerscl(:,(inds(d)-1)*Nhpr+1:(inds(d)-1)*Nhpr+Nhpr)= ...
                            int64(zeros(Ncon*Nhpr,Nhpr));
                    end;
                    
                    % update Rhost
                    Rhostfs=int64(zeros(Ncon,Nhpr));
                    for d=1:Ncon
                        for e=1:Nhpr
                            Rhostfs(d,e)=Rhost(con(d),e)-int64(sum(alloc((d-1)*Nhpr+e,:)));
                        end;
                    end;
                    
                    % re-fairshare with new Rhostfs and R
                    for d=1:Ncon
                        indtemp=find(m{b}(con(d),:)>0);
                        Rhostfs(d,:)=fairshare1(Rhostfs(d,:),R(con(d),indtemp==c)- ...
                            int64(sum(sum(alloc((d-1)*Nhpr+1:(d-1)*Nhpr+Nhpr,:)))));
                        for e=1:Nhpr
                            offerscl((d-1)*Nhpr+e,:)= ...
                                fairshare1(offerscl((d-1)*Nhpr+e,:),Rhostfs(d,e));
                        end;
                    end;
                    
                else
                    chk=0;
                end;
            end;
            
            for d=1:Ncon
                allalloc{con(d)}=[allalloc{con(d)};alloc((d-1)*Nhpr+1:(d-1)*Nhpr+Nhpr,:)];
                allcon{con(d)}=[allcon{con(d)};c];
            end;
            
        end;
        
        
        % ---------- phase 3 ---------- %
        
        % send local traffic indirectly as directed to by receivers
        % alloc has already accounted for rack's sending capacity and
        % host's sending capacity over 1 ToR uplink
        % !!! need to limit sending from each host given multiple pulls
        for c=1:Nrack % sweep racks
            aggD=int64(sum(Dworking((c-1)*Nhpr+1:(c-1)*Nhpr+Nhpr,:,:),3)); % aggregate local demand
            con=find(m{b}(c,:)>0); % get receiver number
            
            % fairshare remaining host capacity
            allalloctemp=int64(zeros(Ncon*Nhpr,Nhost));
            for d=1:Nhpr
                alloctemphost=int64(zeros(Ncon,Nhost));
                for e=1:Ncon
                    alloctemphost(e,:)=allalloc{c}((e-1)*Nhpr+d,:);
                end;
                alloctemphost=fairshare1(alloctemphost(:),Rhost(c,d));
                alloctemphost=reshape(alloctemphost,[Ncon Nhost]);
                for e=1:Ncon
                    allalloctemp((e-1)*Nhpr+d,:)=alloctemphost(e,:);
                end;
            end;
            
            for d=1:Ncon
                ind=find(allcon{c}==con(d));
                alloctemp=allalloctemp((ind-1)*Nhpr+1:(ind-1)*Nhpr+Nhpr,:);
                % we may have been allocated more bytes than we have queued to send
                alloctemp(alloctemp>aggD)=aggD(alloctemp>aggD);
                
                for e=1:Nhpr
                    inds=find(alloctemp(e,:)>0);
                    if isempty(inds)==0 % something was sent
                        [~,nn]=size(inds); % number of elements
                        for f=1:nn
                            [Dsent,Dworking((c-1)*Nhpr+e,inds(f),:),~]= ...
                                fifo(Dworking((c-1)*Nhpr+e,inds(f),:),alloctemp(e,inds(f)));
                            queues{con(d)}((c-1)*Nhpr+e,inds(f),:)= ...
                                queues{con(d)}((c-1)*Nhpr+e,inds(f),:)+Dsent;
                        end;
                    end;
                end;
            end;
        end;
        
    end;
    
%     BWserved(a)=sum(delivered(:))/Nmatch;
    
end;


%%

FCT=(flow_completion-flow_arrival);
FCTcat=FCT(isnan(FCT)==0);

avgFCT=sum(FCTcat)/Nflows;
maxFCT=max(FCTcat);

figure;
plot(delivered(:),'-o','linewidth',2);
xlabel('Matching slot');
ylabel('Total bytes / slot');

%%

save('FCThosts_8racks_16hosts_all2all_10parallel_200kB_Serialized.mat');


















