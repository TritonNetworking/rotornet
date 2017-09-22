
%% Script to plot Figure 9

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


rerunsims=0;
% 1 - rerun the simulations to generate data, then plot
% 0 - load previously saved data, then plot

% NOTE: This simulation needs a lot of RAM if run in parallel
% about 6 GB RAM per worker, so with 12 workers that's ~ 70 GB

%% Run simulation which saves .mat files

path=[pwd '\RotorNet_model_latency\'];

if rerunsims==1
    workers=12; % number of parallel workers, set equal to number of logical processors
    cd RotorNet_model_latency
    Main_RotorLB_256_latency(workers); % this takes a while to run, suggested set workers = 12.
    cd ..
end;

%% load and thin out data

N=256;
Nflows0=(1:200:2400);
Nflows0=[Nflows0 (2400:3000:N^2-N)];
if Nflows0(end)~=N^2-N
    Nflows0=[Nflows0 N^2-N];
end;
NN=max(size(Nflows0)); % number of sweep points
maxsubflow=zeros(1,NN);

for a=1:NN
    fprintf('\nloading file %d of 34\n',a);
    load([path,'subflows',num2str(a),'.mat'],'subflows');
    temp=subflows;
    [~,n]=size(temp);
    if n>0
        maxsubflow(a)=max(temp);
        avgsubflow(a)=mean(temp);
        subflowsthin{a}=temp(1:10000:end);
    else
        subflowsthin{a}=0;
    end;
end;
clear subflows;

%% plot 1

F=18;
LW=2.25;

figure;
hold on;
for a=1:NN
    [~,n]=size(subflowsthin{a});
    h1=plot(Nflows0(a)/Nflows0(end)*ones(1,n),subflowsthin{a},'ok');
end;
grid on;

h2=plot(Nflows0/Nflows0(end),[avgsubflow 0],'-','color','m','linewidth',LW);

h3=plot([0 1],9*[1 1],'--','color','r','linewidth',LW);

hleg=legend([h3 h1 h2],'Delivery bound','Individual subflows','Average');
hleg.Location='NorthEast';

grid on;
box on;
xlim([0 1]);
ylim([0 15]);
set(gca,'fontsize',F);
xlabel('Traffic density');
ylabel('Subflow delivery [slots]');

ax=gca;

ax.XTick=(0:.2:1);

ax.XColor=[0 0 0];
ax.YColor=[0 0 0];
ax.GridColor=[0 0 0];
ax.MinorGridColor=[0 0 0];
ax.GridAlpha=.2;

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%% plot 2

load([path,'subflows30.mat'],'subflows')

figure;
hold on;
h=cdfplot(subflows);
h.LineWidth=LW;
set(gca,'fontsize',F);
xlabel('Subflow delivery time [slots]');
ylabel('CDF');
title('');
box on;
xlim([0 9]);

ax=gca;
ax.YTick=(0:.25:1);


ax.XColor=[0 0 0];
ax.YColor=[0 0 0];
ax.GridColor=[0 0 0];
ax.MinorGridColor=[0 0 0];
ax.GridAlpha=.2;

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width-.05 ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

