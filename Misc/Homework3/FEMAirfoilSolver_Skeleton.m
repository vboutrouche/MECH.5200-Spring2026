close all;
clear all;

%=========================================================================
% FEA/FEM Driver code for the potential flow problem
%=========================================================================




%=========================================================================
% Start by defining your domain: Meshing parameters (done for you)
%=========================================================================
% The following parameters are used for defining the mesh of the object
%
% Shape: This selects an ellispoid or a NACA-4 Digit airfoil
%     -- Ellipsoid: Set the shape to a single value (the aspect ratio of
%     the ellipsoid shape. eg: Shape = [1];
%     -- NACA-4-Digit airfoil: Set the shape parameter a 2-valued entry,
%     eg: Shape = [1 0012], where the second value is the 4-digit reference
%     for the NACA airfoil.
%
% DomainSize: This is a basic parameter that defines how far the domain
% should extend (approx.) from your shape. A value of 4 should be
% sufficient, but you can play around with this.
%
% ref: This is the refinement of the mesh. The higher the value, the more
% elements and vertices you will have.
% 
% powerRef: This is the mesh refinment control. The larger the number, the
% more refined the mesh is near the object. A value of 1 produces a uniform
% refiment. A value of 1.75-2.0 should work well for most of your problems.
%
%=========================================================================

Shape = [1 0045]
DomainSize = 4;
ref = 3;        % <------ Change refinement
powerRef = 2;
[TRI, Nodes, DirichletNodes] = getDiscreteGeometry(Shape, DomainSize, ref, powerRef)

%===
% plot the triangulation
%===
trimesh(TRI, Nodes(:,1),Nodes(:,2))
axis equal
pause(0.1)

%===
% Determine the number of nodes and elements in the domain
%===
NNodes = length(Nodes);
NElem = length(TRI);

%===
% Initialize the A-matrix and RHS to zero valued containers (could use
% spalloc)
%===
A = spalloc(NNodes, NNodes, 12*NNodes);
F = zeros(NNodes, 1);

%===
% cycle through the elements to build the A-matrix and f-vector
%===
for(ie=1:NElem)
    element_number = ie;

    N1 = TRI(ie,1);
    N2 = TRI(ie,2);
    N3 = TRI(ie,3);
    
    X1 = Nodes(N1,1);
    X2 = Nodes(N2,1);
    X3 = Nodes(N3,1);
    
    Y1 = Nodes(N1,2);
    Y2 = Nodes(N2,2);
    Y3 = Nodes(N3,2);
    
    [TriArea] = 1;      % Determine how to calculate and use area

                        % Determine how to compute the planar coefficients
                        
    gN1gN1 = 2;  % Fill the right values
    gN2gN2 = 2;  % Fill the right values
    gN3gN3 = 2;  % Fill the right values

    gN1gN2 = 1;  % Fill the right values
    gN1gN3 = 1;  % Fill the right values             
    gN2gN3 = 1;  % Fill the right values

    A_elemental = TriArea*[gN1gN1  gN1gN2  gN1gN3;...
                   gN1gN2  gN2gN2  gN2gN3;...
                   gN1gN3  gN2gN3  gN3gN3];
               
              
    F_elemental = 0*[1 1 1];    % Fill the RHS
    
    A(N1, N1) = A(N1, N1) + 2; % Fill the right values
    A(N2, N1) = A(N2, N1) + 1; % Fill the right values
    A(N3, N1) = A(N3, N1) + 1; % Fill the right values

    A(N1, N2) = A(N1, N2) + 1; % Fill the right values
    A(N2, N2) = A(N2, N2) + 2; % Fill the right values
    A(N3, N2) = A(N3, N2) + 1; % Fill the right values
    
    A(N1, N3) = A(N1, N3) + 1; % Fill the right values
    A(N2, N3) = A(N2, N3) + 1; % Fill the right values
    A(N3, N3) = A(N3, N3) + 2; % Fill the right values

    F(N1) = F(N1) + 0; % Fill the right values
    F(N2) = F(N2) + 0; % Fill the right values
    F(N3) = F(N3) + 0; % Fill the right values
end

%===
% Boundary conditions
%===
for(i=1:length(DirichletNodes))
    A(DirichletNodes(i),:) = 0;% Fill the right values
    A(DirichletNodes(i),DirichletNodes(i)) = 1;% Fill the right values
    F(DirichletNodes(i)) = 1; % Fill the right values
end

%===
% Solve the problem
%===
Solution = A\F;
Sol = Solution;


%===
% Post Processing calculations: Don't Mess with what works :-)
%===


%===
% Gradient calculation
%===
for(ie=1:NElem)
    element_number = ie;

    N1 = TRI(ie,1);
    N2 = TRI(ie,2);
    N3 = TRI(ie,3);
    
    X1 = Nodes(N1,1);
    X2 = Nodes(N2,1);
    X3 = Nodes(N3,1);
    
    Y1 = Nodes(N1,2);
    Y2 = Nodes(N2,2);
    Y3 = Nodes(N3,2);

    S1 = Sol(N1);
    S2 = Sol(N2);
    S3 = Sol(N3);

    C = [1 X1 Y1; 1 X2 Y2; 1 X3 Y3]\[eye(3)];
    
    Gradient_IE(ie,:) = [(S1*C(2,1)+S2*C(2,2)+S3*C(2,3)), ...
        (S1*C(3,1)+S2*C(3,2)+S3*C(3,3))];
    
    Centroid(ie,:) = [(X1+X2+X3)/3, (Y1+Y2+Y3)/3];
end

%===
% Figure 1
%===
figure
trisurf(TRI, Nodes(:,1), Nodes(:,2), Sol-500) % The 500 is random right now
hold on
trisurf(TRI, Nodes(:,1),-Nodes(:,2),Sol-500)
quiver(Centroid(:,1), Centroid(:,2), Gradient_IE(:,1), Gradient_IE(:,2),.75,'k')
quiver(Centroid(:,1), -Centroid(:,2), Gradient_IE(:,1), -Gradient_IE(:,2),.75,'k')
shading interp
title('Scalar Velocity Potential Distribution (Nodal) and the Velocity Vector Field')
view([0 0 1])
axis equal


%===
% Figure 2
%===
figure
Vel = ((Gradient_IE(:,1).^2 + Gradient_IE(:,2).^2).^(.5))';
a = trisurf(TRI, Nodes(:,1), Nodes(:,2), Nodes(:,2)*0, Vel)
hold on
set(a,'edgealpha',0)
trisurf(TRI, Nodes(:,1), -Nodes(:,2), -Nodes(:,2)*0, Vel)
title('Velocity Distribution (Centroidal)')
view([0 0 1])
axis equal

