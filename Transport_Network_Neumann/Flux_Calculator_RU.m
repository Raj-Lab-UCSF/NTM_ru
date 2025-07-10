%time_qs=0.01;
function [J_0,J_L] = Flux_Calculator_RU(time_qs,varargin)
%%% Solve the  Release-Uptake edge problem with flux boundary conditions 
%%% The outputs are the flux values at x=0 and x=L

close all
format long

%parameters
gamma2_0=0;
frac=0.92;
lambda1_0=1e-02;
lambda2_0=lambda1_0;
delta=1;
epsilon=1e-02;
L1_ = 200;
L2_ = 200;
L_int_ = 1000;
L_ais_ = 40;
L_syn_ = 40;
N1_0_ = 1e-03;
N2_0_ = 5e-03;
M1_0_ = 0;
M2_0_ = 0;

reltol_ = 1e-10;
abstol_ = 1e-10;
fsolvetol_=1e-20;
len_scale_ = 1e-03;
time_scale_ = 1;
va = 0.7*len_scale_ * time_scale_;
vr = 0.7*len_scale_ * time_scale_;
D = 12*len_scale_^2 *time_scale_;
beta = 1e-06*time_scale_;
gamma1_0 = 1e-05 * time_scale_;
gamma1_dt=0e-05*time_scale_;
lambda1_dt=0e-05;

% % New parameters
F_edge_0 = 1e-08 * time_scale_;
mu_r_0 = 1e-05; %release at x=0
mu_r_L = mu_r_0; %release at x=L
mu_u_0 = 1e-04; %uptake at x=0
mu_u_L = mu_u_0; %uptake at x=L
F_edge_dt = 0e-05; %dF/dt
F_edge_dx = 0e-05; %dF/dx


resmesh_ = 'fine';
total_mass_=184;
ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
addParameter(ip, 'beta', beta, validScalar);
addParameter(ip, 'gamma1_0', gamma1_0, validScalar);
addParameter(ip, 'gamma2_0', gamma2_0, validScalar);
addParameter(ip, 'delta', delta, validScalar);
addParameter(ip, 'D', D, validScalar);
addParameter(ip, 'epsilon', epsilon, validScalar);
addParameter(ip, 'F_edge_0',F_edge_0)
addParameter(ip, 'frac', frac, validScalar);
addParameter(ip, 'lambda1_0', lambda1_0, validScalar);
addParameter(ip, 'lambda2_0', lambda2_0, validScalar);
addParameter(ip, 'L_int', L_int_, validScalar);
addParameter(ip, 'L1', L1_, validScalar);
addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'N1_0', N1_0_, validScalar);
addParameter(ip, 'N2_0', N2_0_, validScalar);
addParameter(ip, 'M1_0', M1_0_, validScalar);
addParameter(ip, 'M2_0', M2_0_, validScalar);
addParameter(ip, 'mu_r_0',mu_r_0);
addParameter(ip, 'mu_r_L',mu_r_L);
addParameter(ip, 'mu_u_0',mu_u_0);
addParameter(ip, 'mu_u_L',mu_u_L);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'total_mass',total_mass_);
addParameter(ip, 'va', va);
addParameter(ip, 'vr', vr);
addParameter(ip, 'fsolvetol', fsolvetol_, validScalar);
addParameter(ip, 'reltol', reltol_, validScalar);
addParameter(ip, 'abstol', abstol_, validScalar);
addParameter(ip, 'len_scale', len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);


parse(ip, varargin{:});

L1_new = ip.Results.L1 * ip.Results.len_scale;
L2_new = ip.Results.L2 * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;


L_total = L1_new + L2_new + L_int_new; % size of the system
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000; % number of xmesh points
    num_ext = 200; % number of GM compartments per GM region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,L1_new-10*ip.Results.len_scale,num_ext-40),...
        (L1_new-9.75*ip.Results.len_scale):0.25*ip.Results.len_scale:L1_new];
    xmesh2 = [(L1_new+L_int_new+0.25*ip.Results.len_scale):0.25*ip.Results.len_scale:...
        (L1_new+L_int_new+10*ip.Results.len_scale),...
        linspace(L1_new+L_int_new+10.25*ip.Results.len_scale,L_total,num_ext-40)];
    xmesh_int = [(L1_new+0.25*ip.Results.len_scale):0.25*ip.Results.len_scale:(L1_new+L_ais_new),...
        linspace(L1_new+L_ais_new+0.25*ip.Results.len_scale,L1_new+L_int_new-L_syn_new,...
        num_int-((L_ais_new+L_syn_new)/(0.25*ip.Results.len_scale))),...
        (L1_new+L_int_new-(L_syn_new-0.25*ip.Results.len_scale)):0.25*ip.Results.len_scale:(L1_new+L_int_new)];
elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250; % number of xmesh points
    num_ext = 25; % number of compartments per SD region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,L1_new-10*ip.Results.len_scale,num_ext-5),...
        (L1_new-8*ip.Results.len_scale):2*ip.Results.len_scale:L1_new];
    xmesh2 = [(L1_new+L_int_new+2*ip.Results.len_scale):2*ip.Results.len_scale:...
        (L1_new+L_int_new+10*ip.Results.len_scale),...
        linspace(L1_new+L_int_new+12*ip.Results.len_scale,L_total,num_ext-5)];
    xmesh_int = [(L1_new+2*ip.Results.len_scale):2*ip.Results.len_scale:(L1_new+L_ais_new),...
        linspace(L1_new+L_ais_new+2*ip.Results.len_scale,L1_new+L_int_new-L_syn_new,...
        num_int-((L_ais_new+L_syn_new)/(2*ip.Results.len_scale))),...
        (L1_new+L_int_new-(L_syn_new-2*ip.Results.len_scale)):2*ip.Results.len_scale:(L1_new+L_int_new)];
end
xmesh = [xmesh1, xmesh_int, xmesh2];

presyn_mask2 = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask2);
x1 = xmesh_presyn(end);
ais_mask2 = spatial_mask('ais');
xmesh_ais= xmesh(ais_mask2);
x2 = xmesh_ais(end);
axon_mask2 = spatial_mask('axon');
xmesh_axon= xmesh(axon_mask2);
x3 = xmesh_axon(end);

% % New mesh
xmesh = [xmesh_presyn xmesh_ais xmesh_axon];

gamma1 = @(t) ip.Results.gamma1_0 + t.*gamma1_dt;
lambda1=@(t) ip.Results.lambda1_0 + t.*lambda1_dt;

% % New function 
F_edge = @(x,t) F_edge_0 + t.*F_edge_dt + x.*F_edge_dx; 
%the latter term is only needed to get the vector F_edge_vec in  the
%ode_ss functions 
%replace the term with whatever dependency on x we want

f0 = 0;
t = time_qs;

tic
[~,~,res,flux_0,flux_L]=flux_calculator_RU(xmesh,va,vr,D,t,f0);
time = toc;

J_0 = flux_0;
J_L = flux_L;

fprintf('Residual error %e \n',abs(res))
fprintf('Execution time %e \n', time)


function [n_plot,A_1,res,flux_0,flux_L]=flux_calculator_RU(xmesh,va,vr,D,t,f0)
    B_1=N1_0_; %Extracellular tau in node i
    C_1=N2_0_; %Extracellular tau in node j

    % Shooting parameter: initial value n(x=0)
    % Boundary conditions at the interface x_1 and x_2: continuity of n 

    % Compute the integral of the production term
    int_F_edge = zeros(1,length(xmesh));
    F_edge_vec = F_edge(xmesh,t);
    for j = 2:length(xmesh)
        int_F_edge(j) = trapz(xmesh(1:j),F_edge_vec(1:j));
    end
    int_F_edge = @(x) interp1(xmesh,int_F_edge,x);


    % % % Steady State of Presynaptic Somatodendritic Compartment
    presyn_mask = spatial_mask('presyn');
    xmesh_presyn = xmesh(presyn_mask);
    n0 = @(A) A;
    options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(n0));
    n_ss_presyn = @(A) ode15s(@(x,n)ode_ss_n(x,t,n,A,D),[0,L1_new],n0(A),options);
    n_ss_presyn = @(A,x) deval(n_ss_presyn(A),x);

    % % % Steady state of axon initial segment
    ais_mask = spatial_mask('ais');
    xmesh_ais= xmesh(ais_mask);
    n0 = @(A) n_ss_presyn(A,x1);
    n_ss_ais = @(A) ode15s(@(x,n)ode_ss_n(x,t,n,A,D*lambda1(t)),[L1_new,L1_new+L_ais_new],n0(A),options);
    n_ss_ais = @(A,x) deval(n_ss_ais(A),x);


    % % % Steady state of axon
    axon_mask = spatial_mask('axon');
    xmesh_axon= xmesh(axon_mask);
    n0 = @(A) n_ss_ais(A,x2);
    n_ss_axon = @(A) ode15s(@(x,n)ode_ss_axon(x,t,n,A),[L1_new+L_ais_new,...
        L1_new+L_int_new-L_syn_new],n0(A),options);
    n_ss_axon = @(A,x) deval(n_ss_axon(A),x);
  
    %Calculating the initial data n(0)
    v = @(A,x) ((va*(1+ip.Results.delta.*n_ss_axon(A,x)).*...
        (1-gamma1(t)*ip.Results.epsilon.*n_ss_axon(A,x).^2./(ip.Results.beta-ip.Results.gamma2_0.*n_ss_axon(A,x)))-vr));
    f_init=@(A) -ip.Results.D*ip.Results.frac.*ode_ss_axon(x3,t,n_ss_axon(A,x3),A) + (1-ip.Results.frac).*n_ss_axon(A,x3).*...
       v(A,x3) - mu_r_L*n_ss_axon(A,x3) + mu_u_L*C_1;
    options = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
    A_1=fsolve(f_init,f0,options);
    res = f_init(A_1); %residual error

    n_ss_presyn_plt=n_ss_presyn(A_1,xmesh_presyn);
    n_ss_ais_plt=n_ss_ais(A_1,xmesh_ais);
    n_ss_axon_plt=n_ss_axon(A_1,xmesh_axon);
    n_plot = [n_ss_presyn_plt n_ss_ais_plt n_ss_axon_plt];

    flux_0 = -mu_r_0.*A_1 + mu_u_0.*B_1
    flux_L = mu_r_L.*n_ss_axon(A_1,x3) - mu_u_L.*C_1


    function nprime = ode_ss_n(x,~,~,A,D)
        nprime = 1./D.*(mu_r_0.*A - mu_u_0.*B_1 - int_F_edge(x));
    end

    function nprime = ode_ss_axon(x,t,n,A)
        nprime = 1./(ip.Results.frac*ip.Results.D).*((1-ip.Results.frac).*((va*(1+ip.Results.delta.*n).*...
        (1-((gamma1(t)*ip.Results.epsilon.*n.^2)./(ip.Results.beta-ip.Results.gamma2_0.*n)))-vr)).*n +...
        mu_r_0.*A - mu_u_0.*B_1 - int_F_edge(x));
    end
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

