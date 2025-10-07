function [N, M, R_Edge_Mass] = network_transport_model_Neu_func(matdir,varargin) %(matdir,varargin)

if nargin < 4
    matdir = [cd filesep 'MatFiles'];
 %matdir = [cd filesep 'MatFiles'];
end
alpha_ =0; %1e-2;  % 0; % Recruitment teralm
F_edge_0_ = 2; %1e-8*6*30*24*(60)^2;%1e-08; %Recruitment term edge % 0.5 - 9.9
beta_ = 1e-06;
gamma1_ =0.001; %1e-5 ;%0.001; %2e-03; 
gamma2_ = 0;
delta_ = 50;% 1; %50;
epsilon_ = 25; %1e-02;
    %25;
lambda1_ = 1e-02;   %0.005;%0.025
%lambda2_ = 1e-02; %0.005;% 0.025
gamma1_dt_=0;
lambda_dt_=0;
study_ = 'Hurtado';
init_path_ = [];
init_rescale_ = 2e-2; %2; %2e-2; # Lower rescale
dt_ =0.01;
T_ = 0.1;
%trange_ =[0:0.0005:0.002, 0.0025: 0.0025: 0.1, 0.105:0.005:0.3, 0.31:0.01:2];
%trange_ = [0:0.001:0.004, 0.005:0.005:0.05, 0.06:0.01:0.1, 0.125:0.025:2];
trange_ = [0:0.0001:0.002, 0.003:0.001:0.1, 0.105:0.005:2];

%trange_ = [0:0.0005:0.002, 0.0025:0.0025:0.1, 0.11:0.01:0.2, 0.225:0.025:2];
%trange_ = [0:0.0005:0.002, 0.0025:0.0025:0.1, 0.125:0.025:2];

length(trange_)
%trange_ =[0 0.0005: 0.0025, 0.003: 0.0025: 0.053];
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
%L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'fine';
plotting_ = 1;
reltol_ = 1e-6;
abstol_ = 1e-6;
fsolvetol_ = 1e-20;
connectome_subset_ = 'Hippocampus+PC+RSP';
time_scale_ = 6*30*24*(60)^2;
len_scale_ = 1e-3;

mu_r_0_ = 5; %1e-07*time_scale_; %1.5552 %release at x=0 % 1-9.9
mu_r_L_ = mu_r_0_; %release at x=L
mu_u_0_ = 5; %1e-07*time_scale_; %uptake at x=0 % 1-9.9
mu_u_L_ = mu_u_0_; %uptake at x=L

axon_div_ = 'r1';

ip = inputParser;
% validChar = @(x) ischar(x);
validScalar = @(x) isnumeric(x) && isscalar(x)&& (x>=0);
validLogical = @(x) validScalar(x) && (x == 0 || x == 1);
validAxonDiv=@(x) strcmp(x,'r1') | strcmp(x,'half') | strcmp(x,'r2');
addParameter(ip, 'alpha', alpha_, validScalar);
addParameter(ip, 'F_edge_0', F_edge_0_, validScalar);
addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'gamma1', gamma1_, validScalar);
addParameter(ip, 'gamma2', gamma2_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'lambda1', lambda1_, validScalar);
%addParameter(ip, 'lambda2', lambda2_, validScalar);
addParameter(ip, 'gamma1_dt', gamma1_dt_, validScalar);
addParameter(ip, 'lambda_dt', lambda_dt_, validScalar);
addParameter(ip, 'L_int', L_int_, validScalar);
addParameter(ip, 'L1', L1_, validScalar);
%addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'reltol', reltol_, validScalar);
addParameter(ip, 'abstol', abstol_, validScalar);
addParameter(ip, 'fsolvetol', fsolvetol_, validScalar);
addParameter(ip, 'connectome_subset', connectome_subset_);
addParameter(ip, 'len_scale', len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);
addParameter(ip, 'mu_r_0',mu_r_0_);
addParameter(ip, 'mu_r_L',mu_r_L_);
addParameter(ip, 'mu_u_0',mu_u_0_);
addParameter(ip, 'mu_u_L',mu_u_L_);
addParameter(ip, 'study', study_);
addParameter(ip, 'init_rescale', init_rescale_, validScalar);
addParameter(ip, 'dt', dt_);
addParameter(ip, 'T', T_);
addParameter(ip, 'trange', trange_);
addParameter(ip, 'init_path', init_path_);
addParameter(ip, 'plotting', plotting_, validLogical);
addParameter(ip, 'axon_div', axon_div_, validAxonDiv)
parse(ip, varargin{:});

load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct'); 
load([matdir filesep 'DefaultAtlas.mat'],'DefaultAtlas'); 
load([matdir filesep 'CCF_labels.mat'],'CCF_labels');
load([matdir filesep 'Connectomes.mat'],'Connectomes');
Conn = Connectomes.default;
for i=1:length(Conn)
    for j=1:length(Conn)
       if i==j
           Conn(i,j)=0;
       end
    end
end
thresh_C = 0.8 * mean(nonzeros(Conn(:)));
Conn(Conn < thresh_C) = 0;

Adj = logical(Conn);
% Conn = readmatrix([matdir filesep 'mouse_connectome_19_01.csv']);
% Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);


if ~isempty(ip.Results.init_path)
    init_path = zeros(size(Conn,1),1);
    for i = 1:length(ip.Results.init_path)
        reghemstr = ip.Results.init_path{i};
        reghemcell = split(reghemstr,'_');
        reglog = ismember(CCF_labels(:,1),reghemcell{1});
        if strcmp(reghemcell{2},'L')
            hemlog = ismember(CCF_labels(:,4),'Left Hemisphere'); 
        elseif strcmp(reghemcell{2},'R')
            hemlog = ismember(CCF_labels(:,4),'Right Hemisphere'); 
        else 
            hemlog = ones(size(Conn,1),1);
        end
        init_path((reglog + hemlog) == 2) = 1;
    end
elseif isnan(mousedata_struct.(ip.Results.study).seed)
    %init_path = logical(mousedata_struct.(ip.Results.study).data(:,1));
    init_path = mousedata_struct.(ip.Results.study).data(:,1);
    mean_nz = mean(nonzeros(init_path));
    init_path = init_path ./ mean_nz; % For IbaP301S study where seed is not specified and data is nonbinary at t=1
    init_path = DataToCCF(init_path,ip.Results.study,matdir);
else
    init_path = logical(mousedata_struct.(ip.Results.study).seed);
    init_path = DataToCCF(init_path,ip.Results.study,matdir);
end
beta_new = ip.Results.beta * ip.Results.time_scale;
gamma1_new = ip.Results.gamma1 * ip.Results.time_scale;
gamma2_new = ip.Results.gamma2 * ip.Results.time_scale;
taufun = @(x) ip.Results.init_rescale - (x +...
    (gamma1_new * x.^2)./(beta_new - gamma2_new * x));
options_taufun = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
init_rescale_n =fsolve(taufun,0,options_taufun); %2e-2;
%for i=1:length(init_path)
 id_nan=isnan(init_path);
 init_path(id_nan)=0;
%  end
% endin
%init_path
init_tau = init_rescale_n * init_path;
%init_tau=ip.Results.init_rescale*init_path;
switch ip.Results.connectome_subset
    case 'Hippocampus'
        inds = ismember(CCF_labels(:,3),'Hippocampus');
    case 'Hippocampus+PC+RSP'
        inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
        inds_pc = ismember(CCF_labels(:,1),'Piriform area');
        inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
        inds = logical(inds_hipp + inds_pc + inds_rsp);
    case 'RH'
        inds = ismember(CCF_labels(:,4),'Right Hemisphere');
    case 'LH'
        inds = ismember(CCF_labels(:,4),'Left Hemisphere');
    otherwise
        inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
end
Adj = Adj(inds,inds);
Conn = Conn(inds,inds);
Vol = DefaultAtlas.volumes(inds);
init_tau = init_tau(inds);
nroi = size(Adj,1);
regnamecell = CCF_labels(inds,:);
regnames = cell(size(regnamecell,1),1);
for i = 1:length(regnames)
    regname = regnamecell{i,1};
    reghem = regnamecell{i,4};
    if strcmp(reghem,'Right Hemisphere')
        regnames{i} = [regname ' RH'];
    else
        regnames{i} = [regname ' LH'];
    end
end

if isempty(ip.Results.trange)
    t = 0:ip.Results.dt:ip.Results.T;
else
    t = ip.Results.trange;
end
   nt=length(t);
  N=zeros(nroi,nt);
 
  N(:,1)=init_tau(:);
  % for i=1:nroi
  %     if N(i,1)==0
  %         N(i,1)=1e-4;
  %     end
  % end
    seedregions_ = (N(:,1) > 0);  %    (N(:,1) > 1e-4); % 
 
 i_zero=N(:,1)==0;   %init_tau==0;

Regional_edge_mass = zeros([nroi, nroi, nt, 4]);

 netw_flux_0= zeros([nroi,size(N)]);
  netw_flux_L= zeros([nroi,size(N)]);
  %Mass_edge=zeros([nroi,nroi,nt]);
  n_ss0=zeros([nroi,nroi,nt]);
  %res_max=zeros(1,nt);
 Gamma1=gamma1_new*ones(nroi,nt);
 %Gamma1_der=zeros(nroi,nt);
F_source_edge=zeros([nroi,nroi,nt]);
  gamma_time_dip_roi=[2 6 7 14 2+15 6+15 7+15 14+15];
  %gamma_time_dip_roi=[1:nroi]
  t_gamma=repmat(t, length(gamma_time_dip_roi),1);
%  %size(t_gamma)
 Gamma1(gamma_time_dip_roi, :)=gamma1_new*ones(length(gamma_time_dip_roi),nt) +ip.Results.gamma1_dt*ip.Results.time_scale*t_gamma;
  %Gamma1_der(gamma_time_dip_roi, :)=ip.Results.gamma1_dt*ip.Results.time_scale*ones(length(gamma_time_dip_roi),nt);
Gamma_app=gamma1_new*ones(nroi,1);
 Gamma_app(gamma_time_dip_roi, 1)=gamma1_new*ones(length(gamma_time_dip_roi),1) +ip.Results.gamma1_dt*ip.Results.time_scale*t(1)/2;

  Lambda1=ip.Results.lambda1*ones(nroi,nt);

 lambda_time_dip_roi=[2 6 7 14 2+15 6+15 7+15 14+15];
  %lambda_time_dip_roi=[1:nroi]
  t_lambda=repmat(t, length(lambda_time_dip_roi),1);

 Lambda1(lambda_time_dip_roi, :)=ip.Results.lambda1*ones(length(lambda_time_dip_roi),nt) +ip.Results.lambda_dt*t_lambda;
  %Lambda1_der(lambda_time_dip_roi, :)=ip.Results.lambda_dt*ones(length(lambda_time_dip_roi),nt);
Lambda1_app=ip.Results.lambda1*ones(nroi,1);
 Lambda1_app(gamma_time_dip_roi, 1)=ip.Results.lambda1*ones(length(lambda_time_dip_roi),1) +ip.Results.lambda_dt*t(1)/2;

N_adj_0=[init_rescale_n 0 init_rescale_n]; %N(:,1).*Adj;
N_adj_L=[0 init_rescale_n init_rescale_n];
%Adj
%Idx_netw_x0=logical(N(:,1)).*Adj;
Idx_netw_x0=logical(init_tau).*Adj;

Idx_netw_x0_0= [ 1 1 1];     %logical(N(:,1)).*Adj;
%Idx_netw_x0=logical(Idx_netw_x0+Idx_netw_x0.')
Idx_netw_xL_0=logical([1 1 1].');     %logical(N(:,1));
Idx_netw_xL= logical(init_tau);       %logical(N(:,1));

Gamma1_x0_0=gamma1_new*ones(1,3);     %Gamma1(:,1).*Adj;
%Gamma1_der_x0_0=Gamma1_der(:,1).*Adj;
Lambda1_x0_0=ip.Results.lambda1*ones(1,3);       %Lambda1(:,1).*Adj;
%f_ss0=0.0001*ones(nroi);
f_ss0=[0.0001 0.0001  0.0001];

Adj_in = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);

 fprintf('Flux calculation at initial time \n') 
 tic
 
  [J_0,J_L,~,F_source_edge_0,n_ss_init,~,~]=NetworkFluxCalculator_Neu(N_adj_0,N_adj_L,f_ss0,Adj_in,'beta',ip.Results.beta,...
                                    'delta',ip.Results.delta,'F_edge_0',ip.Results.F_edge_0,...
                                    'epsilon',ip.Results.epsilon,...
                                    'frac',ip.Results.frac,...
                                    'lambda1_x0',Lambda1_x0_0,...
                                    'lambda1_xL',Lambda1_x0_0.',...
                                    'gamma1_x0',Gamma1_x0_0,...
                                    'gamma1_xL',Gamma1_x0_0.',...
                                    'idx_netw_x0',Idx_netw_x0_0,'idx_netw_xL',Idx_netw_xL_0,'frac',ip.Results.frac,...
                                    'L_int',ip.Results.L_int,...
                                    'L1',ip.Results.L1,...
                                    'L_ais',ip.Results.L_ais,...
                                    'L_syn',ip.Results.L_syn,...
                                    'resmesh',ip.Results.resmesh,...
                                    'reltol',ip.Results.reltol,...
                                    'abstol',ip.Results.abstol,...
                                    'fsolvetol',ip.Results.fsolvetol, ...
                                    'axon_div',ip.Results.axon_div,...
                                    'time_scale', ip.Results.time_scale, 'connectome_subset','Single_bis', ...
                                    'mu_r_0',ip.Results.mu_r_0,'mu_r_L',ip.Results.mu_r_L, 'mu_u_0',ip.Results.mu_u_0, 'mu_u_L',ip.Results.mu_u_L ); %compute the steady state  network flux at time t0

  %res_max(1,1);
 netw_flux_0(seedregions_,:,1)=J_0(1,1)*Adj(seedregions_,:);
netw_flux_L(seedregions_,:,1)=J_L(1,1)*Adj(seedregions_,:);
netw_flux_0(:,seedregions_,1)=J_0(1,2)*Adj(:,seedregions_);
netw_flux_L(:,seedregions_,1)=J_L(1,2)*Adj(:,seedregions_);
netw_flux_0(seedregions_,seedregions_,1)=J_0(1,3)*Adj(seedregions_,seedregions_);
netw_flux_L(seedregions_,seedregions_,1)=J_L(1,3)*Adj(seedregions_,seedregions_);
% seedregions_
% Adj(seedregions_,seedregions_)
% netw_flux_0(:,:,1)
% netw_flux_L(:,:,1)
%Mass_edge(seedregions_,:,1)=Mass_edge_0(1,1)*Adj(seedregions_,:);
%Mass_edge(:,seedregions_,1)=Mass_edge_0(1,2)*Adj(:,seedregions_);
%Mass_edge(seedregions_,seedregions_,1)=Mass_edge_0(1,3)*Adj(seedregions_,seedregions_);
F_source_edge(seedregions_,:,1)=F_source_edge_0(1,1)*Adj(seedregions_,:);
F_source_edge(:,seedregions_,1)=F_source_edge_0(1,2)*Adj(:,seedregions_);
F_source_edge(seedregions_,seedregions_,1)=F_source_edge_0(1,3)*Adj(seedregions_,seedregions_);

n_ss0(seedregions_,:,1)=n_ss_init(1,1)*Adj(seedregions_,:);
n_ss0(:,seedregions_,1)=n_ss_init(1,2)*Adj(:,seedregions_);
n_ss0(seedregions_,seedregions_,1)=n_ss_init(1,3)*Adj(seedregions_,seedregions_);

 toc
der=zeros(nroi,nt);

for h=1:  (nt-1)
  
  fprintf('Time step %d/%d\n',h,nt-1) 
 %tic

 %m_t= Gamma1(:,h).*N(:,h).*(2*beta_new-ip.Results.gamma2.*N(:,h)./(beta_new-ip.Results.gamma2.*N(:,h)).^2);
 m_t= Gamma1(:,h).*N(:,h).*((2*beta_new-ip.Results.gamma2.*N(:,h))./(beta_new-ip.Results.gamma2.*N(:,h)).^2);
 
 
 %F_in=time_scale.*netw_flux(:,:,h);
 F_in=netw_flux_L(:,:,h);
 F_out=netw_flux_0(:,:,h);
F_out=F_out.';

 if h<=50 % previously 5
     N_app=N(:,h)+(1./(Vol.*(1+m_t))).*( ( diag((Conn.'*F_in)) - diag((Conn*F_out))) ).*((t(h+1)-t(h))/2);
     %Fun_app=(1./(Vol.*(1+m_t))).*( ( diag((Conn.'*F_in)) - diag((Conn*F_out))) );
     N_adj_app=N_app.*Adj;

  Gamma1_x0_app=Gamma_app.*Adj;
  %Gamma1_der_x0_h1=Gamma1_der(:,h+1).*Adj;
  Lambda1_x0_app=Lambda1_app.*Adj;
[F_out_app,F_in_app,~,~,~,~,~]=NetworkFluxCalculator_Neu(N_adj_app,N_app,n_ss0(:,:,h),Adj_in,'beta',ip.Results.beta,...
                                    'delta',ip.Results.delta,'F_edge_0',ip.Results.F_edge_0,...
                                    'epsilon',ip.Results.epsilon,...
                                    'frac',ip.Results.frac,...
                                    'lambda1_x0',Lambda1_x0_app,...
                                    'lambda1_xL',Lambda1_app,...
                                    'gamma1_x0',Gamma1_x0_app,...
                                    'gamma1_xL',Gamma_app,...
                                    'idx_netw_x0',Idx_netw_x0,'idx_netw_xL',Idx_netw_xL,'frac',ip.Results.frac,...
                                    'L_int',ip.Results.L_int,...
                                    'L1',ip.Results.L1,...
                                    'L_ais',ip.Results.L_ais,...
                                    'L_syn',ip.Results.L_syn,...
                                    'resmesh',ip.Results.resmesh,...
                                    'reltol',ip.Results.reltol,...
                                    'abstol',ip.Results.abstol,...
                                    'fsolvetol',ip.Results.fsolvetol, ...
                                    'axon_div',ip.Results.axon_div,...
                                    'time_scale', ip.Results.time_scale, 'connectome_subset',ip.Results.connectome_subset, ...
                                    'mu_r_0',ip.Results.mu_r_0,'mu_r_L',ip.Results.mu_r_L, 'mu_u_0',ip.Results.mu_u_0, 'mu_u_L',ip.Results.mu_u_L ); %compute the steady state  network flux at time t0
F_out_app=F_out_app.';
%m_t_app=Gamma_app.*N_app.*(2*beta_new-ip.Results.gamma2.*N_app./(beta_new-ip.Results.gamma2.*N_app).^2);
m_t_app=Gamma_app.*N_app.*((2*beta_new-ip.Results.gamma2.*N_app)./(beta_new-ip.Results.gamma2.*N_app).^2);

Fun_app_2=(1./(Vol.*(1+m_t_app))).*( ( diag((Conn.'*F_in_app)) - diag((Conn*F_out_app))) );
N(:,h+1)=N(:,h)+(Fun_app_2).*((t(h+1)-t(h)));
 else
 N(:,h+1)=N(:,h)+(1./(Vol.*(1+m_t))).*( ( diag((Conn.'*F_in)) - diag((Conn*F_out))) ).*(t(h+1)-t(h));  %ip.Results.dt  ; %((diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k);
 end

 N(:,h+1) = N(:,h+1) + ip.Results.alpha*(t(h+1)-t(h))*N(:,h);
 
  N_adj_h1=N(:,h+1).*Adj;
 
  Gamma1_x0_h1=Gamma1(:,h+1).*Adj;
  %Gamma1_der_x0_h1=Gamma1_der(:,h+1).*Adj;
  Lambda1_x0_h1=Lambda1(:,h+1).*Adj;

   fprintf('Flux calculation \n') 

  tic

  % Regional_edge_mass(:,:,h+1,:)
  % Mass_edge(:,:,h+1)

[netw_flux_0(:,:,h+1),netw_flux_L(:,:,h+1),~,F_source_edge(:,:,h+1),n_ss0(:,:,h+1),~,~]=NetworkFluxCalculator_Neu(N_adj_h1,N(:,h+1),n_ss0(:,:,h),Adj_in,'beta',ip.Results.beta,...
                                    'delta',ip.Results.delta,'F_edge_0',ip.Results.F_edge_0,...
                                    'epsilon',ip.Results.epsilon,...
                                    'frac',ip.Results.frac,...
                                    'lambda1_x0',Lambda1_x0_h1,...
                                    'lambda1_xL',Lambda1(:,h+1),...
                                    'gamma1_x0',Gamma1_x0_h1,...
                                    'gamma1_xL',Gamma1(:,h+1),...
                                    'idx_netw_x0',Idx_netw_x0,'idx_netw_xL',Idx_netw_xL,'frac',ip.Results.frac,...
                                    'L_int',ip.Results.L_int,...
                                    'L1',ip.Results.L1,...
                                    'L_ais',ip.Results.L_ais,...
                                    'L_syn',ip.Results.L_syn,...
                                    'resmesh',ip.Results.resmesh,...
                                    'reltol',ip.Results.reltol,...
                                    'abstol',ip.Results.abstol,...
                                    'fsolvetol',ip.Results.fsolvetol, ...
                                    'axon_div',ip.Results.axon_div,...
                                    'time_scale', ip.Results.time_scale, 'connectome_subset',ip.Results.connectome_subset, ...
                                    'mu_r_0',ip.Results.mu_r_0,'mu_r_L',ip.Results.mu_r_L, 'mu_u_0',ip.Results.mu_u_0, 'mu_u_L',ip.Results.mu_u_L ); %compute the steady state  network flux at time t0


 toc

end

 M=(Gamma1.* N.^2)./(beta_new-ip.Results.gamma2 * N);
 R_Edge_Mass = 0; % Regional_edge_mass;
  
%{

[Max_der, I_max_der]=max(der(:,1));

 fprintf('I max der= %d\n',I_max_der)
 M=(Gamma1.* N.^2)./(beta_new-ip.Results.gamma2 * N);
  %M_approx=(Gamma1.* N_approx.^2)./(beta_new-ip.Results.gamma2 * N_approx);

     Mass_node = sum(Vol.*N,1)+sum(Vol.*M,1);
    % size(Mass_node)
     fprintf('Total Node Mass = %d\n',Mass_node)
    Mass_tot_edge=zeros(1,nt);
    for h=1:(nt)
        Mass_tot_edge(1,h) = sum(Conn.*Mass_edge(:,:,h),'all');
    end
    fprintf('Total Edge Mass = %d\n',Mass_tot_edge)
    for h=1:nt
        F=sum(F_source_edge(:,:,h),'all')
    end
     F_source_ed=zeros(1,nt);
    for h=1:nt
      F_source_ed(1,h)=sum(Conn.*F_source_edge(:,:,h),'all').*t(h)
    end
    F_source_nodes=zeros(1,nt);
    source_reg=ones(nroi,1);
    for h=1:nt
    F_source_nodes(1,h)=sum(Vol(seedregions_).*(ip.Results.alpha*source_reg(seedregions_,1)) ,'all').*t(h);
    end
     M_tot=Mass_node+Mass_tot_edge
     F_source=F_source_nodes+F_source_ed
     M_tot_tru=M_tot(1,1)+F_source
     Mass_error=(M_tot-(M_tot(1,1)+F_source))./(M_tot(1,1)+F_source)

    txt = ['$\mathbf{\gamma}_{1,t}$' num2str(ip.Results.gamma1_dt),',', '$\mathbf{\lambda}_{t}$',',' num2str(ip.Results.lambda_dt),',' '$\mathbf{\lambda_1} = $' num2str(lambda1_)  ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
    figure%(1)
    plot(t,M_tot,'r')
    hold on
    plot(t,M_tot_tru,'b')
    ylabel('t')
    xlabel('Mass(t)')
    title('Total mass vs total mass true')
     subtitle(txt,'Interpreter','latex');
%      saveas(figure(1),[ cd '/plot/' 'Tot_mass' '_' 'lambda_var_bis' '.png' ])
%     M_tot(1,1)
%   reldiffs = M_tot - M_tot(1,1);
% reldiffs = reldiffs ./ M_tot(1,1);
% 
figure %(2); 
%hold on;
%for i = 1:size(masstots,1)
    plot(t,Mass_error); 
    ylabel('t')
    xlabel('E(t)')
    title('Relative error Total mass')
 subtitle(txt,'Interpreter','latex');
%     saveas(figure(2),[ cd '/plot/' 'Rel_err' '_' 'lambda_var_bis' '.png' ])
% %end
% 
% 
% 
% 
figure%(3)
%figure
subplot(2,1,1)
plot(t,N);


xlabel('t');
ylabel('N(t)');
 title("N,M distributions on the network",'Fontsize',12);
  %txt = ['$\mathbf{\gamma}_{1,t} = $' num2str(gamma1_dt), ',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 %subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M);
xlabel('t');
ylabel('M(t)');
%  saveas(figure(3),[ cd '/plot/' 'N_M' '_' 'lambda_var_bis' '.png' ])
% 

figure%(4)
%figure
subplot(2,1,1)
plot(t,N(i_zero,:));
 title("N,M distributions on the network",'Fontsize',12);
 % txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',', '$\mathbf{\lambda}_{t}$',',' num2str(lambda_dt),',' '$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 %subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M(i_zero,:));
%  saveas(figure(4),[ cd '/plot/' 'N_M_no_seed' '_' 'lambda_var_bis' '.png' ])

figure%(5)
%figure

plot(t,res_max);


xlabel('t');
ylabel('res_max');
 title("max residual shooting parameter calculation",'Fontsize',12);


% 
% figure(5)
% %figure
% plot(t,N+M);
% 
% 
% xlabel('t');
% ylabel('N(t)+M(t)');
%  title("N,M distributions on the network",'Fontsize',12);
%   %txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
%  subtitle(txt,'Interpreter','latex');
%  saveas(figure(5),[ cd '/plot/' 'Tot_tau' '_' 'lambda_var_bis' '.png' ])
% 
%  %figure(5)
% figure(6)
% plot(t,N(i_zero,:)+M(i_zero,:));
% 
% 
% xlabel('t');
% ylabel('N(t)+M(t)');
%  title("Total Tau on the network",'Fontsize',12);
%   %txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
%  subtitle(txt,'Interpreter','latex');
%  saveas(figure(6),[ cd '/plot/' 'Tot_tau_no_seed' '_' 'lambda_var_bis' '.png' ])
% %nroi_idx=1:nroi;
%   %figure(6)
% figure(7)
% heatmap(N(i_zero,1:end)+M(i_zero,1:end));
% 
% 
% xlabel('t');
% ylabel('N(t)+M(t)');
%  title("Total Tau on the network");
% % subtitle(txt,'Interpreter','latex');
%  saveas(figure(7),[ cd '/plot/' 'Tot_tau_heatmap' '_' 'lambda_var' '.png' ])
%  %  txt = ['$\mathbf{\gamma}_{1,t}$' num2str(gamma1_dt),',','$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\beta}=$' num2str(beta_) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
%  % subtitle(txt,'Interpreter','latex');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% figure(8)
% for h=1:nt
% 
%     subplot(2,1,1)
%     %heatmap(1:nroi,1:nroi, B(:,:,h))
%     plot(t(h),max(abs(B(:,:,h)),[],'all'),'o',LineWidth=2)
%     hold on
%     title('B network')
%     %subtitle(txt,'Interpreter','latex')
%     subplot(2,1,2)
%     %heatmap(1:nroi,1:nroi,6*30*24*(60)^2*netw_flux(:,:,h))
% 
%     plot(t(h),max(abs(6*30*24*(60)^2*netw_flux(:,:,h)),[],'all'),'*',LineWidth=2)
%    hold on
%  % xlabel('nroi')
%  % ylabel('J,B (o-* line)' )
%  title('J network')
% 
% end
% % saveas(figure(8),[ cd '/plot/' 'B_J' '_' 'lambda_var' '.png' ])
% saveas(figure(8),[ cd '/plot/' 'BJ' '_' 'lambda_var_bis' '.png' ])

end

%}