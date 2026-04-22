%[text] # Week 11 Activity: Transient Heat Equation — Quadratic Elements, Forward Euler
%[text] Build a 1D FEM solver for the transient heat equation using **quadratic (second-order)** elements and **forward Euler** (explicit) time integration with a **lumped mass** matrix. Complete each **TODO** block.
%%
%[text] ## 0) Problem Statement
%[text] We solve the transient heat equation:
%[text]{"align":"center"} $ \frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2}, \quad x \in (0,1),\; t > 0 $
%[text] with homogeneous Dirichlet boundary conditions $ u(0,t) = u(1,t) = 0 $ and initial condition $ u(x,0) = \sin(\pi x) $.
%[text] The exact solution is $ u(x,t) = e^{-\pi^2 t}\sin(\pi x) $.

L       = 1.0;
t_final = 0.1;

u_exact = @(x,t) exp(-pi^2*t) .* sin(pi*x);
u_init  = @(x) sin(pi*x);

x_plot = linspace(0, L, 300);
figure
plot(x_plot, u_init(x_plot), 'LineWidth', 1.8);
grid on; xlabel('x'); ylabel('u(x,0)');
title('Initial Condition: u(x,0) = sin(\pi x)');
%%
%[text] ## 1) Semi-Discrete Weak Form
%[text] Multiply by test function $ v $, integrate by parts:
%[text]{"align":"center"} $ \int_0^1 \frac{\partial u}{\partial t}\,v\,dx + \int_0^1 \frac{\partial u}{\partial x}\frac{\partial v}{\partial x}\,dx = 0 $
%[text] This yields the semi-discrete system $ M\,\dot{U} + K\,U = 0 $.
%%
%[text] ## 2) Quadratic Basis Functions
%[text] Three-node reference element on $ \xi \in [-1,1] $.

syms xi_s real

%%% TODO 1 — Define the three quadratic shape functions
%  Hint: N1 = xi*(xi-1)/2,  N2 = 1-xi^2,  N3 = xi*(xi+1)/2
N1s = ____;
N2s = ____;
N3s = ____;
%%% END TODO 1

N = [N1s; N2s; N3s];

xi_num = linspace(-1, 1, 200);
figure
plot(xi_num, double(subs(N1s, xi_s, xi_num)), 'LineWidth', 1.8); hold on;
plot(xi_num, double(subs(N2s, xi_s, xi_num)), 'LineWidth', 1.8);
plot(xi_num, double(subs(N3s, xi_s, xi_num)), 'LineWidth', 1.8);
grid on; xlabel('\xi'); ylabel('Value');
title('Quadratic Basis Functions');
legend('N_1 (left)', 'N_2 (mid)', 'N_3 (right)', 'Location', 'bestoutside');
%%
%[text] ## 3) Element Matrices — Symbolic Derivation
%[text] **Stiffness (3x3):** $ k_{ab}^{(e)} = \int_{-1}^{1} \frac{dN_a}{dx}\frac{dN_b}{dx}\,J\,d\xi $
%[text] **Lumped mass:** $ m_a^{(e)} = \int_{-1}^{1} N_a\,J\,d\xi $

syms h_e positive
J    = h_e/2;
dNdx = (1/J) * diff(N, xi_s);

%%% TODO 2 — Compute element stiffness symbolically (3x3)
ke_sym = sym(zeros(3));
for a = 1:3
    for b = 1:3
        ke_sym(a,b) = ____;
    end
end
ke_sym = simplify(ke_sym)
%%% END TODO 2

%%% TODO 3 — Compute lumped mass vector symbolically (3x1)
me_lumped_sym = sym(zeros(3,1));
for a = 1:3
    me_lumped_sym(a) = ____;
end
me_lumped_sym = simplify(me_lumped_sym)
%%% END TODO 3
%%
%[text] ## 4) Assembly
%[text] With $ n_{el} $ quadratic elements the total node count is $ 2 n_{el}+1 $. Connectivity: [left, mid, right].

nel     = 10;
nn      = 2*nel + 1;
x_nodes = linspace(0, L, nn).';
conn    = [(1:2:2*nel-1).', (2:2:2*nel).', (3:2:2*nel+1).'];

K      = zeros(nn);
M_diag = zeros(nn, 1);

for e = 1:nel
    nd = conn(e,:);
    he = x_nodes(nd(3)) - x_nodes(nd(1));

    %%% TODO 4 — Numeric element stiffness and lumped mass
    ke   = ____;   % 3x3
    me_l = ____;   % 3x1
    %%% END TODO 4

    %%% TODO 5 — Scatter into global matrices
    K(nd, nd)  = ____;
    M_diag(nd) = ____;
    %%% END TODO 5
end

figure; spy(K);
xlabel('j'); ylabel('i');
title('Stiffness K — sparsity (quadratic)');
%%
%[text] ## 5) Initial Condition

U = u_init(x_nodes);
U(1) = 0; U(end) = 0;
%%
%[text] ## 6) Forward Euler Time Integration
%[text]{"align":"center"} $ U^{n+1} = U^n - \Delta t\,M_L^{-1}\,K\,U^n $
%[text] **Stability:** for quadratic elements with lumped mass, $ \Delta t < h_e^2/7 $ (approximately).

he = L / nel;
dt = 0.8 * he^2 / 7;
Nt = ceil(t_final / dt);
dt = t_final / Nt;
fprintf('dt = %.6f,  Nt = %d\n', dt, Nt);

t_snap = [0, t_final/2, t_final];
U_snap = zeros(nn, numel(t_snap));
U_snap(:,1) = U;
snap_idx = 2;

for n = 1:Nt
    tn = n * dt;

    %%% TODO 6 — Forward Euler update
    %  Hint: R = -K*U;  U = U + dt*(R ./ M_diag);  then enforce BCs
    R = ____;
    U = ____;
    U(1) = 0;  U(end) = 0;
    %%% END TODO 6

    if snap_idx <= numel(t_snap) && tn >= t_snap(snap_idx) - dt/2
        U_snap(:,snap_idx) = U;
        snap_idx = snap_idx + 1;
    end
end
%%
%[text] ## 7) Validation

figure; hold on;
colors = lines(numel(t_snap));
leg = {};
for k = 1:numel(t_snap)
    plot(x_nodes, U_snap(:,k), 'o-', 'Color', colors(k,:), ...
         'LineWidth', 1.5, 'MarkerSize', 4);
    fplot(@(x) u_exact(x, t_snap(k)), [0 L], '--', ...
         'Color', colors(k,:), 'LineWidth', 1.8);
    leg{end+1} = sprintf('FEM  t=%.3f', t_snap(k));
    leg{end+1} = sprintf('Exact t=%.3f', t_snap(k));
end
grid on; xlabel('x'); ylabel('u(x,t)');
title('Forward Euler — Quadratic Elements');
legend(leg, 'Location', 'best');

err = max(abs(U - u_exact(x_nodes, t_final)));
fprintf('Max nodal error at t=%.2f: %.4e\n', t_final, err);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
