function [network_flux_0, network_flux_L, mass_edge,F_source_edge,n_ss0,res_max,region_edge_mass] = NetworkFluxCalculator_Neu(tau_x0,tau_xL,f_ss0,Adj_input,varargin) % removed matdir

%if nargin < 3
%    matdir = [cd filesep 'MatFiles'];
    %matdir=['/Users/veronicatora/Desktop/ODEJustin/'];
%end

% % % 1. Preset values of flexible parameters
beta_ = 1e-06;
gamma1_x0_= 2e-05;
gamma1_xL_ = 2e-05 ;
gamma2_ = 0;
delta_ = 50;%      
epsilon_ =25;    
lambda1_x0_ = 0.0025;    
lambda1_xL_ = 0.0025; 
% lambda1_der_x0_ = 0;%0.01; %0.01;%0.02  0.01 
% lambda1_der_xL_ =0 ; %0.01;
% gamma1_der_x0_ =0; % 0.01; %0.01;%0.02  0.01 
% gamma1_der_xL_ =0; 
%lambda2_=0.025;
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
%L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'coarse';
reltol_ = 1e-6;
abstol_ = 1e-6;
fsolvetol_ = 1e-6;
connectome_subset_ ='Hippocampus+PC+RSP'; %'Hippocampus';
len_scale_ = 1e-3;
time_scale_ = 1; %6*30*24*(60)^2% 1;
% % New parameters
F_edge_0_ = 1e-08 * time_scale_;
mu_r_0_ = 1e-05* time_scale_; %release at x=0
mu_r_L_ = mu_r_0_; %release at x=L
mu_u_0_ = 1e-04* time_scale_; %uptake at x=0
mu_u_L_ = mu_u_0_; %uptake at x=L
F_edge_dt = 0e-05; %dF/dt
F_edge_dx = 0e-05; %dF/dx
idx_netw_x0_=1;
idx_netw_xL_=1;
axon_div_='r1';



ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validArray = @(x) isnumeric(x); % && (sum(x>=0) == length(x));
validLogical=@(x) islogical(x);
validAxonDiv=@(x) strcmp(x,'r1') | strcmp(x,'half') | strcmp(x,'r2');

addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'F_edge_0', F_edge_0_ , validScalar);
addParameter(ip, 'gamma2', gamma2_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'gamma1_x0', gamma1_x0_, validArray);
addParameter(ip, 'gamma1_xL', gamma1_xL_, validArray);
addParameter(ip, 'lambda1_x0', lambda1_x0_, validArray);
addParameter(ip, 'lambda1_xL', lambda1_xL_, validArray);
 addParameter(ip, 'idx_netw_x0', idx_netw_x0_, validArray);
 addParameter(ip, 'idx_netw_xL', idx_netw_xL_, validLogical);
% addParameter(ip, 'gamma1_der_x0', gamma1_der_x0_, validArray);
% addParameter(ip, 'gamma1_der_xL', gamma1_der_xL_, validArray);

addParameter(ip, 'mu_r_0',mu_r_0_,validScalar);
addParameter(ip, 'mu_r_L',mu_r_L_,validScalar);
addParameter(ip, 'mu_u_0',mu_u_0_,validScalar);
addParameter(ip, 'mu_u_L',mu_u_L_,validScalar);
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
addParameter(ip, 'axon_div', axon_div_, validAxonDiv)
parse(ip, varargin{:});
beta_new = ip.Results.beta*ip.Results.time_scale;
%gamma1_new = ip.Results.gamma1*ip.Results.time_scale;
%gamma2_new = ip.Results.gamma2*ip.Results.time_scale;
L1_new = ip.Results.L1 * ip.Results.len_scale;
%L2_new = ip.Results.L2 * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;
gamma1_x0 =ip.Results.gamma1_x0 ;
gamma1_xL =ip.Results.gamma1_xL;
lambda1_x0 =ip.Results.lambda1_x0 ;
lambda1_xL =ip.Results.lambda1_xL;
idx_netw_x0= ip.Results.idx_netw_x0;
idx_netw_xL=ip.Results.idx_netw_xL;

% gamma1_der_x0 =ip.Results.gamma1_der_x0 ;
% gamma1_der_xL =ip.Results.gamma1_der_xL;
% lambda1_der_x0 =ip.Results.lambda1_der_x0 ;
% lambda1_der_xL =ip.Results.lambda1_der_xL;

% % % 2. Definition of static constants
v_a = 0.7*ip.Results.len_scale * ip.Results.time_scale; % Average velocity (um/s) of anterograde active transpot (Konsack 2007)
v_r = 0.7*ip.Results.len_scale * ip.Results.time_scale; % Average velocity (um/s) of retrograde active transport (Konsack 2007)
diff_n = 12*ip.Results.len_scale^2 * ip.Results.time_scale; % Diffusivity (um^2/s) of n (Konsack 2007)
mu_r_0 = ip.Results.mu_r_0; %release at x=0
mu_r_L = ip.Results.mu_r_L; %release at x=L
mu_u_0 = ip.Results.mu_u_0; %uptake at x=0
mu_u_L = ip.Results.mu_u_L; %uptake at x=L


% % % 3. Definition of the (inhomogeneous) xmesh
L_total = L1_new + L_int_new -L_syn_new; % size of the system
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000; % number of xmesh points
    num_ext = 100; % number of GM compartments per GM region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,L1_new-10*ip.Results.len_scale,num_ext-40),...
        (L1_new-9.75*ip.Results.len_scale):0.25*ip.Results.len_scale:L1_new];
    % xmesh2 = [(L1_new+L_int_new+0.25*ip.Results.len_scale):0.25*ip.Results.len_scale:...
    %     (L1_new+L_int_new+10*ip.Results.len_scale),...
    %     linspace(L1_new+L_int_new+10.25*ip.Results.len_scale,L_total,num_ext-40)];
    xmesh_int = [(L1_new+0.25*ip.Results.len_scale):0.25*ip.Results.len_scale:(L1_new+L_ais_new),...
        linspace(L1_new+L_ais_new+0.25*ip.Results.len_scale,L1_new+L_int_new-L_syn_new,...
        num_int-((L_ais_new+L_syn_new)/(0.25*ip.Results.len_scale))),...
        (L1_new+L_int_new-(L_syn_new-0.25*ip.Results.len_scale)):0.25*ip.Results.len_scale:(L1_new+L_int_new-L_syn_new)];
elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250; % number of xmesh points
    num_ext = 25; % number of compartments per SD region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,L1_new-10*ip.Results.len_scale,num_ext-5),...
        (L1_new-8*ip.Results.len_scale):2*ip.Results.len_scale:L1_new];
    % xmesh2 = [(L1_new+L_int_new+2*ip.Results.len_scale):2*ip.Results.len_scale:...
    %     (L1_new+L_int_new+10*ip.Results.len_scale),...
    %     linspace(L1_new+L_int_new+12*ip.Results.len_scale,L_total,num_ext-5)];
    xmesh_int = [(L1_new+2*ip.Results.len_scale):2*ip.Results.len_scale:(L1_new+L_ais_new),...
        linspace(L1_new+L_ais_new+2*ip.Results.len_scale,L1_new+L_int_new-L_syn_new,...
        num_int-((L_ais_new+L_syn_new)/(2*ip.Results.len_scale))),...
        (L1_new+L_int_new-(L_syn_new-2*ip.Results.len_scale)):2*ip.Results.len_scale:(L1_new+L_int_new-L_syn_new)];
end
xmesh = [xmesh1, xmesh_int];



% % % 4. Steady State Calculation
% Follows the derivation of Michiel Bertsh
%gamma_1 time dipendent convex combination
% gamma1_fun=@(gamma1_i,gamma1_j,x)(1-x./L_total).*(gamma1_i)+(x./L_total).*(gamma1_j);
% gamma1_der_fun=@(gamma1_der_i,gamma1_der_j,x)(1-x./L_total).*(gamma1_der_i)+(x./L_total).*(gamma1_der_j);
% %lambda_fun=@(lambda_i,lambda_j,x)(1-x./L_total).*(lambda_i)+(x./L_total).*(lambda_j);
gamma1_fun=@(gamma1_i,gamma1_j,x)(gamma1_i+gamma1_j)./2;
lambda_fun=@(lambda_i,lambda_j,x)(lambda_i+lambda_j)./2;
% F_edge = @(x,t) (F_edge_0 + t.*F_edge_dt + x.*F_edge_dx);
%  int_F_edge = zeros(1,length(xmesh));
%  F_edge_vec = F_edge(xmesh,t);
%  for j = 2:length(xmesh)
%         int_F_edge(j) = trapz(xmesh(1:j),F_edge_vec(1:j));
%  end
%  int_F_edge = @(x,idx) interp1(xmesh,int_F_edge,x).*idx;
 x0=xmesh(1);
  int_F_edge = @(x,idx) (ip.Results.F_edge_0.*(x-x0)).*idx;
% % % 4a. Presynaptic somatodendritic compartment
presyn_mask = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask);
n0 = @(B) B;
options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol);    %,'NonNegative'  ,1:length(n0));
n_ss_presyn = @(B,N_i,idx) ode15s(@(x,n)ode_ss_n(x,B,n,diff_n,N_i,idx),[0,L1_new],n0(B),options);
n_ss_presyn = @(B,N_i,idx,x) deval(n_ss_presyn(B,N_i,idx),x);
x1 = xmesh_presyn(end);

% % % 4b. Axon initial segment
ais_mask = spatial_mask('ais');
xmesh_ais = xmesh(ais_mask);
x2 = xmesh_ais(end);
n_ss_ais = @(B,N_i,idx,lambda1_i,lambda1_j) ode15s(@(x,n)ode_ss_n(x,B,n,diff_n*lambda_fun(lambda1_i,lambda1_j,x),N_i,idx),[L1_new,L1_new+L_ais_new],n_ss_presyn(B,N_i,idx,x1),options);
n_ss_ais = @(B,N_i,idx,lambda1_i,lambda1_j,x) deval(n_ss_ais(B,N_i,idx,lambda1_i,lambda1_j),x);


% % % 4c. Axon
axon_mask = spatial_mask('axon');
xmesh_axon= xmesh(axon_mask);
n_ss_axon = @(B,N_i,idx,lambda1_i,lambda1_j,gamma1_i,gamma1_j) ode15s(@(x,n)ode_ss_axon(x,B,gamma1_i,gamma1_j,n,N_i,idx),[L1_new+L_ais_new,...
   L1_new+L_int_new-L_syn_new],n_ss_ais(B,N_i,idx,lambda1_i,lambda1_j,x2),options);
n_ss_axon = @(B,N_i,idx,lambda1_i,lambda1_j,gamma1_i,gamma1_j,x) deval(n_ss_axon(B,N_i,idx,lambda1_i,lambda1_j,gamma1_i,gamma1_j),x);
x3 = xmesh_axon(end);
 %xmesh(end)
 f_init=@(B,N_i,N_j,idx,lambda1_i,lambda1_j,gamma1_i,gamma1_j) -mu_r_0.*B + mu_u_0.*N_i+ int_F_edge(x3,idx)- mu_r_L.*n_ss_axon(B,N_i,idx,lambda1_i,lambda1_j,gamma1_i,gamma1_j,x3) + mu_u_L.*N_j;
 % % % 5a. Flux calculation on network 
 %ones(length(tau_x0))

 %Adj=readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
Adj = Adj_input;

switch ip.Results.connectome_subset
    case 'Hippocampus'
        Adj = Adj([27:37 (27+213):(37+213)], [27:37 (27+213):(37+213)]);
    case 'Hippocampus+PC+RSP'
        adjinds = [27:37,78:80,147];
        adjinds = [adjinds,adjinds+213];
        Adj = Adj(adjinds,adjinds);
    case 'RH'
        Adj = Adj(1:213,1:213);
    case 'LH'
        Adj = Adj(214:end,214:end);
    case 'Single'
        Adj = 1;
    case 'Single_bis'
        Adj=[1 1 1]; %if there are connections between two different seedregions otherwise Adj=[1 1]
    
end
nroi = size(Adj,2);

network_flux_0 = zeros(size(Adj));
network_flux_L = zeros(size(Adj));

%Gamma1_i=zeros(size(Adj));
%Gamma1_j=zeros(size(Adj));

%Lambda1_i=0.01*ones(size(Adj));
%Lambda1_j=0.01*ones(size(Adj));
n_ss0=zeros(size(Adj));
res=zeros(size(Adj));
F_source_edge=zeros(size(Adj));

n_ss0_dict = zeros(size(Adj));
res_dict = zeros(size(Adj));
F_source_edge_dict = zeros(size(Adj));
network_flux_0_dict = zeros(size(Adj));
network_flux_L_dict = zeros(size(Adj));

parfor i = 1:nroi

    Adj_in = logical(Adj(:,i));
    tau_xL_i = tau_xL(i);
    tau_xL_i = repmat(tau_xL_i,length(Adj_in),1);
     idx_netw_xL_i=idx_netw_xL(i);
     idx_netw_xL_i = repmat(idx_netw_xL_i,length(Adj_in),1);
    gamma1_xL_i = gamma1_xL(i);
    gamma1_xL_i = repmat(gamma1_xL_i,length(Adj_in),1);
    
     lambda1_xL_i = lambda1_xL(i);
    lambda1_xL_i = repmat(lambda1_xL_i,length(Adj_in),1);

    
    i_app =Adj_in;
    tau_x0_i = tau_x0(i_app,i);
    tau_xL_i = tau_xL_i(i_app);
    idx_netw_x0_i= idx_netw_x0(i_app);
    idx_netw_xL_i= idx_netw_xL_i(i_app);
    gamma1_x0_i = gamma1_x0(i_app,i);
    gamma1_xL_i = gamma1_xL_i(i_app);
    %Gamma1_i(i_app,i)=gamma1_x0_i;
    %Gamma1_j(i_app,i)=gamma1_xL_i;
    
    lambda1_x0_i = lambda1_x0(i_app,i);
    lambda1_xL_i = lambda1_xL_i(i_app);
    %Lambda1_i(i_app,i)=lambda1_x0_i;
    %Lambda1_j(i_app,i)=lambda1_xL_i;
    
    if ~isempty(tau_x0_i)

        options = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
        idx_netw=logical(idx_netw_x0_i+idx_netw_xL_i);
        fun_ss = @(B) f_init(B,tau_x0_i,tau_xL_i,idx_netw,lambda1_x0_i,lambda1_xL_i,gamma1_x0_i,gamma1_xL_i);   
       
        f0=f_ss0(i_app,i)+1e-7;

        %n_ss0(i_app,i)=fsolve(fun_ss,f0,options);

        n_ss0_temp = fsolve(fun_ss,f0,options);
        
        temp_n_ss0_arr = zeros(nroi,1);
        temp_n_ss0_arr(i_app) = n_ss0_temp;
        n_ss0_dict(i,:)=temp_n_ss0_arr;

        %res(i_app,i)=fun_ss(n_ss0(i_app,i));

        res_temp = fun_ss(n_ss0_temp);
        temp_res_arr = zeros(nroi,1);
        temp_res_arr(i_app) = res_temp;
        res_dict(i,:) = temp_res_arr;

        %F_source_edge(i_app,i)= int_F_edge(x3,idx_netw);

        F_source_edge_temp = int_F_edge(x3,idx_netw);
        temp_F_source_edge = zeros(nroi,1);
        temp_F_source_edge(i_app) = F_source_edge_temp;
        F_source_edge_dict(i,:)= temp_F_source_edge;

        %network_flux_0(i_app,i) = -mu_r_0.* n_ss0(i_app,i)+ mu_u_0.*tau_x0_i;

        netflux_0_temp = -mu_r_0.* n_ss0_temp+ mu_u_0.*tau_x0_i;
        temp_netflux_0_arr = zeros(nroi,1);
        temp_netflux_0_arr(i_app) = netflux_0_temp;
        network_flux_0_dict(i,:) = temp_netflux_0_arr;

        %network_flux_L(i_app,i)=mu_r_L.*n_ss_axon(n_ss0(i_app,i),tau_x0_i,idx_netw,lambda1_x0_i,lambda1_xL_i,gamma1_x0_i,gamma1_xL_i,x3) - mu_u_L.*tau_xL_i;
        
        netflux_L_temp = mu_r_L.*n_ss_axon(n_ss0_temp,tau_x0_i,idx_netw,lambda1_x0_i,lambda1_xL_i,gamma1_x0_i,gamma1_xL_i,x3) - mu_u_L.*tau_xL_i;
        temp_netflux_L_arr = zeros(nroi,1);
        temp_netflux_L_arr(i_app) = netflux_L_temp;
        network_flux_L_dict(i,:) = temp_netflux_L_arr;

    end
 end

for i = 1:nroi

    Adj_in = logical(Adj(:,i));
    i_app =Adj_in;

    n_ss0(i_app,i) = n_ss0_dict(i,i_app);
    res(i_app,i) = res_dict(i,i_app);
    F_source_edge(i_app,i) = F_source_edge_dict(i,i_app);
    network_flux_0(i_app,i) = network_flux_0_dict(i,i_app);
    network_flux_L(i_app,i) = network_flux_L_dict(i,i_app);

end

% %%% Calculate edge mass belonging to region 1 and region 2

mass_edge = 0;
region_edge_mass = 0;

% axon_div = ip.Results.axon_div;
% 
% presyn_len = size(xmesh_presyn,2);
% ais_len = size(xmesh_ais,2);
% 
% if strcmp(axon_div,'r2')
% 
%     n_ss_r1 = [n_ss_presyn_eval n_ss_ais_eval];
%     n_ss_r2 = [n_ss_axon_eval];
% 
%     xmesh_r1 = [xmesh_presyn xmesh_ais];
%     xmesh_r2 = [xmesh_axon];
% 
% elseif strcmp(axon_div,'half')
% 
%     mid_index = fix(size(xmesh_axon,2) / 2);
% 
%     n_ss_r1 = [n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval(:,1:mid_index)];
%     n_ss_r2 = [n_ss_axon_eval(:,mid_index:end)];
% 
%     xmesh_r1 = [xmesh_presyn xmesh_ais xmesh_axon(:,1:mid_index)];
%     xmesh_r2 = [xmesh_axon(:,mid_index:end)];
% 
% else % strcmmp(axon_div,'r1') - default
% 
%     n_ss_r1 = [n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval];
%     n_ss_r2 = [];
% 
%     xmesh_r1 = [xmesh_presyn xmesh_ais xmesh_axon];
%     xmesh_r2 = [];
% end
% 
% m_ss_r1 = (gamma1_fun(Gamma1_i,Gamma1_j,[xmesh_r1]).* n_ss_r1.^2)./(beta_new-ip.Results.gamma2 *n_ss_r1);
% 
% r1_n_mass = trapz(xmesh_r1, n_ss_r1, 2);
% r1_m_mass = trapz(xmesh_r1, m_ss_r1, 2);
% 
% if strcmp(axon_div,'r1')
%     r2_n_mass = zeros(size(r1_n_mass));
%     r2_m_mass = zeros(size(r1_m_mass));
% else
%     m_ss_r2 = (gamma1_fun(Gamma1_i,Gamma1_j,[xmesh_r2]).* n_ss_r2.^2)./(beta_new-ip.Results.gamma2 *n_ss_r2);
%     r2_n_mass = trapz(xmesh_r2, n_ss_r2, 2);
%     r2_m_mass = trapz(xmesh_r2, m_ss_r2, 2);
% end
% 
% r1_n_mass = reshape(r1_n_mass,size(Adj));
% r1_m_mass = reshape(r1_m_mass,size(Adj));
% r2_n_mass = reshape(r2_n_mass,size(Adj));
% r2_m_mass = reshape(r2_m_mass,size(Adj));
% 
% region_edge_mass(:,:,1) = r1_n_mass;
% region_edge_mass(:,:,2) = r1_m_mass;
% region_edge_mass(:,:,3) = r2_n_mass;
% region_edge_mass(:,:,4) = r2_m_mass;

% figure
% plot(max(abs(res),[],1))

res_max=max(abs(res),[],'all');
%res_max=0;


% % % 6. Functions
   
    function nprime = ode_ss_n(x,B,n,D,N_i,idx)
        nprime = (1./D).*(mu_r_0.*B - mu_u_0.*N_i - int_F_edge(x,idx));
    end

    function nprime = ode_ss_axon(x,B,gamma1_i,gamma1_j,n,N_i,idx)
        
        nprime = 1./(ip.Results.frac*diff_n).*((1-ip.Results.frac).*((v_a*(1+ip.Results.delta.*n).*...
        (1-((gamma1_fun(gamma1_i,gamma1_j)*ip.Results.epsilon.*n.^2)./(beta_new-ip.Results.gamma2.*n)))-v_r)).*n +...
        mu_r_0.*B - mu_u_0.*N_i - int_F_edge(x,idx));
    end

    function [maskvals] = spatial_mask(compartment)
        switch compartment
            case 'presyn'
                maskvals = (xmesh <= L1_new);
            case 'ais'
                maskvals = logical(-1 + (xmesh > L1_new) + ...
                    (xmesh < (L1_new + L_ais_new)));
            case 'axon'
                maskvals = logical(-1 + (xmesh >= L1_new + L_ais_new) + ...
                    (xmesh <= (L1_new + L_int_new - L_syn_new)));
            case 'syncleft'
                maskvals = logical(-1 + (xmesh > L1_new + L_int_new - L_syn_new) + ...
                    (xmesh < (L1_new + L_int_new))); 
            case 'postsyn'
                maskvals = logical(-1 + (xmesh >= L1_new + L_int_new) + ...
                    (xmesh <= (L1_new + L_int_new + L2_new)));
            otherwise
                error('Incorrect compartment specification')
        end
    end

end
