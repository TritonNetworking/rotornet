function []=Main_RotorLB_testbed()

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

N=8; % # network endpoints
Nreal=32; % number of demand matrix realizations of each # flows
Ncycles=4;
filename='RotorLB_testbed_32real.mat'; % .mat file to save output

%% set up matchings

Nmatch=4; % number of matchings per switch (crossover matchings)
NS=ceil(N/Nmatch); % number of Rotor switches

% the bidirectional matchings used in the testbed:
mc=cell(1,log2(N));
for a=1:log2(N)
    mc{a}=zeros(N);
    for b=1:N/(2^a)
        for c=0:(2^a-1)
            mc{a}((2^a)*b-c,(2^a)*b-((2^a-1)-c))=1;
        end;
    end;
end;

mc{4}=mc{3};
mc{3}=zeros(N);
mc{3}(1,6)=1;
mc{3}(2,5)=1;
mc{3}(3,8)=1;
mc{3}(4,7)=1;
mc{3}(5,2)=1;
mc{3}(6,1)=1;
mc{3}(7,4)=1;
mc{3}(8,3)=1;

m=cell(1,Nmatch); % all connections at each of Nmatch timeslots
for a=1:Nmatch
    m{a}=mc{a};
    m{a}=m{a}+fliplr(eye(N))*mc{a};
    % ensure main diagonal is zero
    m{a}(logical(diag(ones(1,N),0)))=zeros(N,1);
end;

%% Main loop:

path=[pwd '\Demand_matrices\'];

Nflows0=(1:1:N^2-N); % setup to sweep number of non-zero flows
NN=max(size(Nflows0)); % number of sweep points

workers=8; % number of parallel workers
parpool(workers); % start parallel pool

[BWr_1hop,BWr_2hop]=deal(zeros(NN,Nreal));

parfor x=1:Nreal
    for y=1:NN
        
        fprintf(sprintf('\n\nIteration: y = %d, x = %d',y,x));
        
        file=sprintf('Nelems=%d_Realization=%d.txt',y,x);
        D=csvread([path file]);
        
        % ----------- get modeled RotorNet solutions ----------- %
        
        % only 1-hop paths
        BWserved=RotorLB_1hop(N,Nmatch,NS,m,1e6*D,Ncycles); % get BW from RotorLB
        BWr_1hop(y,x)=BWserved(4);
        
        % 1-hop AND 2-hop paths
        BWserved=RotorLB_2hop(N,Nmatch,NS,m,1e6*D,Ncycles); % get BW from RotorLB
        BWr_2hop(y,x)=BWserved(4);
        
    end;
end;

delete(gcp); % close parallel workers

% save variables:
save(filename);

%%







