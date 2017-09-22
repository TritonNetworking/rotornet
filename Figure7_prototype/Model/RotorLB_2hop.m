function [BWserved]=RotorLB_2hop(N,Nmatch,NS,m,D,Ncycles)

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

%% initialize queues
queues=cell(N,3); % (endpoint number, [original,indirect])
for a=1:N
    queues{a,1}=D(a,:); % direct queues
    queues{a,2}=zeros(N); % indirect queues
    queues{a,3}=zeros(N); % temp queues
end;

%% Simulator:

tol=1e-8; % tolerance to round extra link capacity

BWserved=zeros(1,Ncycles);

for a=1:Ncycles % iterate through cycles of matchings
    
    delivered=zeros(N); % initialize delivered bandwidth for this cycle
    
    for b=1:Nmatch % iterate through matchings
        
        % initialize temp queues - these will be used to update indirect
        % queues
        for x=1:N
            queues{x,3}=zeros(N);
        end;
        
        % initialize link capacities
        con=find(m{b}(1,:)>0);
        [~,Ncon]=size(con); % get number of connections
        linkcap=1/NS; % fraction of ToR BW in each RotorNet link
        R=linkcap*ones(N,Ncon);
        
        % ---------- phase 1 ---------- %
        
        % first, serve any existing nonlocal (two hop) traffic
        for c=1:N % sweep local nodes
            con=find(m{b}(c,:)>0); % get ids of connected nodes
            for d=1:Ncon % sweep receivers
                vserved=fairshare1(queues{c,2}(:,con(d)),R(c,d)); % unused link
                
                % update delivered bandwidth
                delivered(:,con(d))=delivered(:,con(d))+vserved;
                % update nonlocal queues
                queues{c,2}(:,con(d))=queues{c,2}(:,con(d))-vserved;
                % remaining capacity on the link
                R(c,d)=R(c,d)-sum(vserved);
                
            end;
        end;
        
        % next, send local data directly if possible
        for c=1:N % sweep local nodes
            con=find(m{b}(c,:)>0); % get ids of connected nodes
            for d=1:Ncon
                if R(c,d)>0 % there is remaining capacity on the link
                    if R(c,d)>queues{c,1}(con(d)) % there is more than enough capacity
                        delivered(c,con(d))=delivered(c,con(d))+queues{c,1}(con(d)); % update delivered
                        R(c,d)=R(c,d)-queues{c,1}(con(d)); % update remaining capacity
                        queues{c,1}(con(d))=0; % update local queues
                    else
                        delivered(c,con(d))=delivered(c,con(d))+R(c,d); % update delivered
                        queues{c,1}(con(d))=queues{c,1}(con(d))-R(c,d); % update local queues
                        R(c,d)=0; % no more capacity on the link
                    end;
                end;
            end;
        end;
        
        allalloc=cell(1,N); % store allocation arrays by sending node
        % each cell has Ncon entries
        allcon=cell(1,N); % store which node is allocating
        % each cell has Ncon entries
        
        % bidirectional step - each endpoint computes how much indirect
        % data it can accept from each sender, given its own queues and the senders' offers
        for c=1:N % sweep local nodes
            con=find(m{b}(c,:)>0); % get ids of connected nodes
            % get local queues of each connected node - giving how much traffic the sender want to indirect
            % fairshare that according to how much link capacity remains
            offer=zeros(Ncon,N);
            for d=1:Ncon
                offer(d,:)=queues{con(d),1};
            end;
            bufferavail=linkcap*ones(1,N)-(queues{c,1}+sum(queues{c,2}));
            bufferavail(bufferavail<0)=0;
            
            inds=find(bufferavail==0); % remove any offers where buffer = 0
            [~,Ninds]=size(inds);
            for d=1:Ninds
                offer(:,inds(d))=zeros(Ncon,1);
            end;
            
            offerscl=zeros(Ncon,N);
            
            for d=1:Ncon
                indtemp=find(m{b}(con(d),:)>0);
                offerscl(d,:)=fairshare1(offer(d,:),R(con(d),indtemp==c));
            end;
            
            alloc=zeros(Ncon,N); % bandwidth allocated to senders
            chk=1;
            while chk==1
                if Ncon>1
                    slices=find(sum(offerscl)>0); % which indirect queues to consider fairsharing over
                else
                    slices=find(offerscl>0);
                end;
                [~,Nslices]=size(slices);
                if Nslices>0
                    for d=1:Nslices
                        temp=fairshare1(offerscl(:,slices(d)),bufferavail(slices(d)));
                        bufferavail(slices(d))=round((bufferavail(slices(d))-sum(temp))/tol)*tol;
                        alloc(:,slices(d))=alloc(:,slices(d))+temp;
                    end;
                    offerscl=offer-alloc;
                    
                    inds=find(bufferavail<=0); % remove any offers where buffer = 0
                    [~,Ninds]=size(inds);
                    for d=1:Ninds
                        offerscl(:,inds(d))=zeros(Ncon,1);
                    end;
                    
                    for d=1:Ncon
                        indtemp=find(m{b}(con(d),:)>0);
                        offerscl(d,:)=fairshare1(offerscl(d,:),R(con(d),indtemp==c)-sum(alloc(d,:)));
                    end;
                else
                    chk=0;
                end;
            end;
            
            for d=1:Ncon
                allalloc{con(d)}=[allalloc{con(d)};alloc(d,:)];
                allcon{con(d)}=[allcon{con(d)};c];
            end;
            
        end;
        
        
        % ---------- phase 2 ---------- %
        
        % simply send local traffic indirectly as directed to by receivers
        for c=1:N % sweep local nodes
            con=find(m{b}(c,:)>0); % get receiver number
            for d=1:Ncon
                alloctemp=allalloc{c}(allcon{c}==con(d),:);
                % we may have been allocated more BW than we have queued to send
                alloctemp(alloctemp>queues{c,1})=queues{c,1}(alloctemp>queues{c,1});
                
                % update local queue
                queues{c,1}=queues{c,1}-alloctemp;
                
                % update temp queues
                queues{con(d),3}(c,:)=queues{con(d),3}(c,:)+alloctemp;
            end;
        end;
        
        % update receiving indirect queues from temp queues
        for c=1:N
            queues{c,2}=queues{c,2}+queues{c,3};
        end;
        
    end;
    
    BWserved(a)=sum(delivered(:))/Nmatch;
    
end;

%%












