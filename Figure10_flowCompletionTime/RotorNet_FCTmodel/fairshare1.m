function vfs=fairshare1(v0,alloc)

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

v=v0;

velem=v>0;
Nelem=sum(velem);

if ~Nelem==0
    chk=0;
    while chk==0
        f=idivide(alloc,int64(Nelem),'floor'); % fairshare allocation
        v=v-f*int64(velem);
        if min(v)<0 % some elements were overserved
            elemtemp=v<0;
            alloc=-1*int64(sum(v(elemtemp)));
            v(elemtemp)=0;
            velem=v>0;
            Nelem=sum(velem);
            if Nelem==0
                chk=1;
            end;
        else % no elements were overserved
            chk=1;
        end;
    end;
end;

vfs=v0-v;

% dif=alloc-int64(sum(vfs)) % number of unused bytes
% if dif>0
%     inds=randi([1 Nelem],1,dif);
%     for a=1:dif
%         vfs
%     end;
% end;