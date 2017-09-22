function []=Main_RotorLB_16()

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

N=16; % # network endpoints
Nreal=93; % number of demand matrix realizations of each # flows
Ncycles=4;
filename='RotorLB_N=16_93.mat'; % .mat file to save output

%% set up matchings

% assign number of switches (NS) and number of matchings per switch
% (Nmatch):

Nmatch=4; % number of matchings per switch (crossover matchings)
NS=ceil(N/Nmatch); % number of switches

% construct bidirectional matchings
mc=cell(1,log2(N));
for a=1:log2(N)
    mc{a}=zeros(N);
    for b=1:N/(2^a)
        for c=0:(2^a-1)
            mc{a}((2^a)*b-c,(2^a)*b-((2^a-1)-c))=1;
        end;
    end;
end;

C{1}=[1 2 3];
C{2}=[1 2 3 4];
C{3}=[4];

[~,NC]=size(C);

m=cell(1,Nmatch); % all connections seen at each of Nmatch timeslots
for a=1:Nmatch
    m{a}=mc{a};
    for b=1:NC
        [~,n]=size(C{b});
        temp=eye(N);
        for c=1:n
            temp=temp*mc{C{b}(c)};
        end;
        m{a}=m{a}+temp*mc{a};
    end;
    % ensure main diagonal is zero
    m{a}(logical(diag(ones(1,N),0)))=zeros(N,1);
end;

%% Main loop:

Nflows0=(1:10:N^2-N); % setup to sweep number of non-zero flows
if Nflows0(end)~=N^2-N
    Nflows0=[Nflows0 N^2-N];
end;

NN=max(size(Nflows0)); % number of sweep points

workers=8; % number of parallel workers
parpool(workers); % start parallel pool

BWr=zeros(NN,Nreal,Ncycles);

parfor y=1:Nreal
    for x=1:NN
        
        fprintf(['\n\nIteration: y = ',num2str(y),', x = ',num2str(x)]);
        
        % generate a random demand matrix:
        
        r=[ones(1,Nflows0(x)) zeros(1,N^2-N-Nflows0(x))]; % fill array
        r=r(randperm(length(r))); % randomly shuffle
        % set diag to zeros
        D=zeros(N); % initialize
        for b=0:N-1
            s=N*b-b+1; % start
            D(b+1,:)=[r(s:s+b-1) 0 r(s+b:s+N-2)];
        end;
        
        D=1e3*D; % scale so we don't run out of demand
        
        % ----------- get RotorNet solutions ----------- %
        
        BWr(x,y,:)=RotorLB(N,Nmatch,NS,m,D,Ncycles); % get bandwidth delivered by RotorLB
        
    end;
end;

delete(gcp); % close parallel workers

% save variables:
save(filename);


%%



