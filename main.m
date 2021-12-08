clear ; tic; clc
warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Codes for Example 1 in paper:
%%% Localized Collocation Schemes: Various Formulations and Recent Applications
%%% Copyright: Zhuojia Fu's group in Hohai University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameters
dim = 3;
AN = 10;
n = AN^3;
N = 6*AN^2;
ns = 50;

%% computational domain
dx = 1/(AN+1);
ttt = dx:dx:1-dx;
[BX,BY] = meshgrid(ttt,ttt);
BX = reshape(BX,1,AN^2);
BY = reshape(BY,1,AN^2);
A(:,1:6*AN^2) = [ones(1,AN^2) zeros(1,AN^2) BY         BX          BX          BY;...
    BX      BY        ones(1,AN^2) zeros(1,AN^2)  BY          BX;...
    BY      BX        BX         BY          ones(1,AN^2)   zeros(1,AN^2)];
[X,Y,Z] = meshgrid(ttt,ttt,ttt);
A(:,6*AN^2+1:n+6*AN^2) = [reshape(X,1,n);reshape(Y,1,n);reshape(Z,1,n)];
nou1=[ones(1,AN^2) -ones(1,AN^2) zeros(1,AN^2) zeros(1,AN^2) zeros(1,AN^2) zeros(1,AN^2);.../
    zeros(1,AN^2) zeros(1,AN^2) ones(1,AN^2) -ones(1,AN^2) zeros(1,AN^2) zeros(1,AN^2);.../
    zeros(1,AN^2) zeros(1,AN^2) zeros(1,AN^2) zeros(1,AN^2) ones(1,AN^2) -ones(1,AN^2)];
f = A.';
nou=nou1.';

%% the vector of b in Cx=b
tn=6*N/6;
exu=f(:,1).^2-f(:,2).^2+f(:,3);
uf1(1:tn,1)=f(1:tn,1).^2-f(1:tn,2).^2+f(1:tn,3);
uf1(N+1:N+n,1)=0;

%% LCM (Localized Collocation Methods)
%%
[E,ix] = GFDM_3D(f,ns,N+n); % GFDM
%%
% [E,ix] = LRBFCM_3D(f,ns,N+n,5); % LRBFCM
%%
% [E,~,~,~,ix] = LMFS_3D(f,ns,N+n,dim); % LMFS and LMFS1
% E(1,1,:) = E(1,1,:)-1;
%%
% [E,~,~,~,ix] = LRTCM_3D(f,ns,N+n,dim,0.5); % LRTCM and LRTCM1
% E(1,1,:) = E(1,1,:)-1;
%%
% [E,ix] = LCTM_3D(f,ns,N+n); % LCTM and LCTM1
% E(1,1,:) = E(1,1,:)-1;

%% discrete sparse matrix C
C = sparse(N+n,N+n); 
for i=1:tn
    C(i,i) = 1;
end

for i=tn+1:N+n
    for j=1:ns
        a = ix(j,i);
        C(i,a) = E(1,j,i);
    end
end

uu = C\uf1; 

Rerr = norm(abs(uu(1:N+n,1)-exu(1:N+n,1)))/norm(exu(1:N+n,1));
disp(['  Total number of points   ','    ：',num2str(N+n)]);
disp(['  Global error   ','：',num2str(Rerr)]);
disp(['  CPU time   ','  ：',num2str(toc),' s']);
beep,pause(0.5),beep,pause(0.5),beep;
