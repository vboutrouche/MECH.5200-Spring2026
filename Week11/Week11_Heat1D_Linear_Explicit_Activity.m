%[text] # Week 11 Activity: Transient Heat Equation — Linear Elements, Forward Euler
%[text] Build a 1D FEM solver for the transient heat equation using **linear (first-order)** elements and **forward Euler** (explicit) time integration with a **lumped mass** matrix. Complete each **TODO** block.
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
%[text] Substituting $ u_h = \sum_j U_j(t)\,N_j(x) $ and choosing $ v = N_i $ yields the **semi-discrete system**:
%[text]{"align":"center"} $ M\,\dot{U} + K\,U = 0 $
%[text] where $ M_{ij} = \int N_i\,N_j\,dx $ (mass matrix) and $ K_{ij} = \int N_i'\,N_j'\,dx $ (stiffness matrix).
%%
%[text] ## 2) Linear Basis Functions
%[text] Two-node reference element on $ \xi \in [-1,1] $.

syms xi_s real

%%% TODO 1 — Define the two linear shape functions N1s, N2s
%  Hint: N_1(xi) = (1-xi)/2,  N_2(xi) = (1+xi)/2
N1s = ____;
N2s = ____;
%%% END TODO 1

N = [N1s; N2s];

xi_num = linspace(-1, 1, 200);
figure
plot(xi_num, double(subs(N1s, xi_s, xi_num)), 'LineWidth', 1.8); hold on;
plot(xi_num, double(subs(N2s, xi_s, xi_num)), 'LineWidth', 1.8);
grid on; xlabel('\xi'); ylabel('Value');
title('Linear Basis Functions');
legend('N_1', 'N_2', 'Location', 'best');
%%
%[text] ## 3) Element Matrices — Symbolic Derivation
%[text] Jacobian $ J = h_e/2 $, derivatives $ dN_a/dx = (1/J)\,dN_a/d\xi $.
%[text] **Stiffness:** $ k_{ab}^{(e)} = \int_{-1}^{1} \frac{dN_a}{dx}\frac{dN_b}{dx}\,J\,d\xi $
%[text] **Lumped mass:** $ m_a^{(e)} = \int_{-1}^{1} N_a\,J\,d\xi $

syms h_e positive
J    = h_e/2;
dNdx = (1/J) * diff(N, xi_s);

%%% TODO 2 — Compute the 2x2 element stiffness matrix symbolically
ke_sym = sym(zeros(2));
for a = 1:2
    for b = 1:2
        ke_sym(a,b) = ____;   % Hint: integrate dNdx(a)*dNdx(b)*J over xi_s in [-1,1]
    end
end
ke_sym = simplify(ke_sym)
%%% END TODO 2

%%% TODO 3 — Compute the lumped mass vector symbolically
me_lumped_sym = sym(zeros(2,1));
for a = 1:2
    me_lumped_sym(a) = ____;  % Hint: integrate N(a)*J over xi_s in [-1,1]
end
me_lumped_sym = simplify(me_lumped_sym)
%%% END TODO 3
%%
%[text] ## 4) Assembly
%[text] Assemble global stiffness $ K $ and lumped mass vector $ M_d $.

NumElems = 20;
NumNodes = NumElems + 1;
xNodes   = linspace(0, L, NumNodes).';

K      = zeros(NumNodes);
M_diag = zeros(NumNodes, 1);

for e = 1:NumElems
    nodes_e = [e, e+1];
    he = xNodes(e+1) - xNodes(e);

    %%% TODO 4 — Write the numeric element stiffness and lumped mass
    ke   = ____;   % 2x2 matrix — substitute h_e = he into ke_sym
    me_l = ____;   % 2x1 vector — substitute h_e = he into me_lumped_sym
    %%% END TODO 4

    %%% TODO 5 — Scatter into global matrices
    K(nodes_e, nodes_e) = ____;
    M_diag(nodes_e)     = ____;
    %%% END TODO 5
end

figure; spy(K);
xlabel('j'); ylabel('i');
title('Stiffness K — sparsity');
%%
%[text] ## 5) Initial Condition

U = u_init(xNodes);
U(1) = 0; U(end) = 0;
%%
%[text] ## 6) Forward Euler Time Integration
%[text] Using lumped mass the update is fully explicit:
%[text]{"align":"center"} $ U^{n+1} = U^n - \Delta t\,M_L^{-1}\,K\,U^n $
%[text] **Stability (CFL):** for uniform linear elements, $ \Delta t < h^2/2 $.

he = L / NumElems;
dt = 0.8 * he^2 / 2;           % 80 % of CFL limit
Nt = ceil(t_final / dt);
dt = t_final / Nt;
fprintf('dt = %.6f,  Nt = %d,  CFL ratio = %.2f\n', dt, Nt, dt/(he^2/2));

t_snap = [0, t_final/2, t_final];
U_snap = zeros(NumNodes, numel(t_snap));
U_snap(:,1) = U;
snap_idx = 2;

for n = 1:Nt
    tn = n * dt;

    %%% TODO 6 — Forward Euler update
    %  Hint: R = -K*U;  U = U + dt*(R ./ M_diag);
    %  Then enforce BCs: U(1) = 0; U(end) = 0;
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
    plot(xNodes, U_snap(:,k), 'o-', 'Color', colors(k,:), ...
         'LineWidth', 1.5, 'MarkerSize', 4);
    fplot(@(x) u_exact(x, t_snap(k)), [0 L], '--', ...
         'Color', colors(k,:), 'LineWidth', 1.8);
    leg{end+1} = sprintf('FEM  t=%.3f', t_snap(k));
    leg{end+1} = sprintf('Exact t=%.3f', t_snap(k));
end
grid on; xlabel('x'); ylabel('u(x,t)');
title('Forward Euler — Linear Elements');
legend(leg, 'Location', 'best');

err = max(abs(U - u_exact(xNodes, t_final)));
fprintf('Max nodal error at t=%.2f: %.4e\n', t_final, err);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
