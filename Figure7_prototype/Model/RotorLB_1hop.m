function [BWserved]=RotorLB_1hop(N,Nmatch,NS,m,D,Ncycles)

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
        
    end;
    
    BWserved(a)=sum(delivered(:))/Nmatch;
    
end;

%%












