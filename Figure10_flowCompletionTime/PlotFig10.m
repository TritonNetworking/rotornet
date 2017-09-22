
%% Script to plot Figure 10

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


%% load data from ns3 simulation

% 16 hosts, 8 racks, all-to-all inter-rack (240 flows) x 10 flows parallel

fid=fopen('FCT_FatTree.csv','r');
data=textscan(fid,'%s','delimiter',',','headerlines',1);
data=data{1};

strt=18817; % Fat tree, 240 flows, 10 parallel flows, 200 kB / flow
stp=23296;

cnt=0;
for a=strt:stp
    dstport=str2double(sscanf(data{(a-1)*21+21},'%s'));
    if dstport==8080
        cnt=cnt+1;
        fctft(cnt)=str2double(sscanf(data{(a-1)*21+7},'%s'))- ...
            str2double(sscanf(data{(a-1)*21+4},'%s'));
    end;
end;

fctft=fctft*1e3; % convert to ms


%% Run RotorNet simulation and load data

path=[pwd '\RotorNet_FCTmodel\'];

% Run RotorNet simulation function
if rerunsims==1
    cd RotorNet_FCTmodel
    RotorLB_8racks_RLhost_FCT_FS(); % fair share BW between flows
    RotorLB_8racks_RLhost_FCT_Serial(); % serialize BW to flows
    cd ..
end;


%% plot

F=16; % fontsize
LW=2; % linewidth


figure;
hold on;
set(gcf,'position',[572 784 792 420]);

%--------- Fat Tree

h1=cdfplot(fctft);
h1.LineWidth=LW;
h1.Color='r';

y=get(h1,'ydata');
x=get(h1,'xdata');

plot(x(1:1000:end),y(1:1000:end),'xr','markersize',10,'linewidth',LW);
h1=plot(-1,-1,'-xr','linewidth',LW,'markersize',10);

ylim([0 1]);
xlim([0 45]);

val=.999; % 99.9th percentile
ind=find(y>=val,1);
plot(x(ind)*[1 1],[0 1],'--','linewidth',LW,'color','r');

xlims=get(gca,'XLim');
ylims=get(gca,'YLim');
x1frac=.89;
y1frac=.24;
ht = text((1-x1frac)*xlims(1)+x1frac*xlims(2), ...
    (1-y1frac)*ylims(1)+y1frac*ylims(2),'99.9%-tile Fat Tree');
set(ht,'Rotation',90);
set(ht,'FontSize',F);
ht.Color='r';

%--------- RotorNet fair

load([path 'FCThosts_8racks_16hosts_all2all_10parallel_200kB.mat']);

FCTcat=1e3*FCTcat; % ms

h2=cdfplot(FCTcat);
h2.LineWidth=LW;
h2.Color=.5*[0 1 0];

y=get(h2,'ydata');
x=get(h2,'xdata');

plot(x(1:150:end),y(1:150:end),'o','markersize',10,'linewidth',LW,'color',.5*[0 1 0]);
h2=plot(-1,-1,'-or','linewidth',LW,'markersize',10,'color',.5*[0 1 0]);

%--------- RotorNet unfair

load([path 'FCThosts_8racks_16hosts_all2all_10parallel_200kB_Serialized.mat']);

FCTcat=1e3*FCTcat; % ms

h3=cdfplot(FCTcat);
h3.LineWidth=LW;
h3.Color=.6*[0 1 1];

y=get(h3,'ydata');
x=get(h3,'xdata');

plot(x(1:50:end),y(1:50:end),'d','markersize',10,'linewidth',LW,'color',.6*[0 1 1]);
h3=plot(-1,-1,'-dr','linewidth',LW,'markersize',10,'color',.6*[0 1 1]);

val=1; % 100th percentile
ind=find(y>=val,1);
plot(x(ind)*[1 1],[0 1],'--k','linewidth',LW);

xlims=get(gca,'XLim');
ylims=get(gca,'YLim');
x1frac=.71;
y1frac=.2;
ht = text((1-x1frac)*xlims(1)+x1frac*xlims(2), ...
    (1-y1frac)*ylims(1)+y1frac*ylims(2),{'100%-tile SelectorNet';'     (fair & unfair)'});
set(ht,'Rotation',90);
set(ht,'FontSize',F);
ht.Color='k';

hleg=legend([h1 h2 h3],'Fat Tree','SelectorNet (fair)','SelectorNet (unfair)');
hleg.Location='northwest';
hleg.FontSize=F;

ax=gca;
ax.FontSize=F;
ylabel('CDF');
xlabel('Flow completion time (ms)');
title('');
grid on;
box on;

ax.XColor=[0 0 0];
ax.GridColor=[0 0 0];
ax.MinorGridColor=[0 0 0];
ax.GridAlpha=.2;

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width-.01 ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


%%



