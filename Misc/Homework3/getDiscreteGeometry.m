
function [TRI, Nodes, DirichletNodes] = getDiscreteGeometry(Shape, DomainSize, ref, powerRef)


%===
% This code is a Willis hack to define the geometry for a symmetric 
% potential flow solver that usesFEM
%===
XX = DomainSize; 
YY = DomainSize;
PowerRefine = powerRef;


OuterBoundary = [-XX, 0; -XX, YY; 1/2 YY; 1+XX YY; 1+XX 0];

%========================================================================
% Refine the outer boundary 
%========================================================================
for(j=1:(ref-1))
NewBoundary = [];
    for(i=1:(length(OuterBoundary(:,1))-1))
        % Refine outer boundary -- cycle through boundary elements to
        % divide them 
        NewBoundary = [NewBoundary; OuterBoundary(i,:)];
        NewBoundary = [NewBoundary; (OuterBoundary(i,:)+OuterBoundary(i+1,:))/2];
    end
    NewBoundary = [NewBoundary; OuterBoundary(end,:)];    
    OuterBoundary = NewBoundary;
end

if(length(Shape)==1)
    CircleNodes = (length(OuterBoundary)-1);
    theta = 0:1:(CircleNodes);
    theta = theta/max(theta)*pi;
    InnerBoundary = [-cos(theta)/2+1/2; sin(theta)/2/Shape]';
    InnerBoundary(end,2) = InnerBoundary(end,2)+1e-10;
else
    CircleNodes = (length(OuterBoundary)-1);
    theta = 0:1:(CircleNodes);
    theta = (theta.^1.1)/max(theta.^1.1)*pi;
    xval = -cos(theta)/2+1/2;
    
    params.plot_on = 0;
	NACA = Shape(2);
	[ x_up, y_up, x_lo, y_lo, x_cl, y_cl ] = getNACA(NACA, xval, params);
    InnerBoundary = [x_up; y_up]'
    InnerBoundary(end,2) = InnerBoundary(end,2)+1e-10;

end

%========================================================================
% Set nodes between the two boundaries
%========================================================================
Nodes = [];
for(j = 1:(ref*4+1))
    tNodes = [];
    for(i = 1:length(InnerBoundary(:,1)))
       tNodes = [tNodes; ((j-1)/(ref*4))^PowerRefine*(OuterBoundary(i,:) - InnerBoundary(i,:)) + InnerBoundary(i,:)];       
    end
    Nodes = [Nodes; tNodes];
    % Storing the boundaries in the desired format
    if(j==1)
        Boundary = Nodes;
    elseif(j>1 && j<(ref*4+1))        
        Boundary = [Boundary; [tNodes(1,:); tNodes(end,:)]];
    else
        Boundary = [Boundary; OuterBoundary];
    end
end

%========================================================================
% Define Closed Boundary of domain
%========================================================================
InnerBoundaryRev = [];
for(i=1:length(InnerBoundary(:,1)))
    InnerBoundaryRev(i,:) = InnerBoundary(end-i+1,:);
end
Outline = [OuterBoundary; InnerBoundaryRev; OuterBoundary(1,:)];

%========================================================================
% Use Delaunary to triangulate the domain
%========================================================================
TRI = delaunay(Nodes);

%========================================================================
% Cycle through triangles and eliminate exterior to shape triangles
%========================================================================
TRINOW = [];
for(i=1:length(TRI))
    centroid = [sum(Nodes(TRI(i,1:3),1))/3, sum(Nodes(TRI(i,1:3),2))/3];
    TotalAngle = 0;
    for(j=2:length(Outline))
       TempAngle1 = centroid - Outline(j-1,:);
       TempAngle2 = centroid - Outline(j,:);
       TempAngle1(3) = 0;
       TempAngle2(3) = 0;       
       Dangle = cross(TempAngle1, TempAngle2)./norm(TempAngle1)./norm(TempAngle2);
       TotalAngle = TotalAngle + asin(Dangle(3));       
    end    
    if((abs(TotalAngle))<(pi-eps))
    else
        TRINOW = [TRINOW; TRI(i,:)];
    end
end
TRI = TRINOW;

TRINOW = [];
for(i=1:length(TRI(:,1)))
    
    N1 = TRI(i,1);
    N2 = TRI(i,2);
    N3 = TRI(i,3);
    
    X1 = Nodes(N1,1);
    X2 = Nodes(N2,1);
    X3 = Nodes(N3,1);
    
    Y1 = Nodes(N1,2);
    Y2 = Nodes(N2,2);
    Y3 = Nodes(N3,2);
    
    [TriArea] = calculateArea(X1, X2, X3, Y1, Y2, Y3);
    
    if(TriArea>1e-15)
        TRINOW = [TRINOW; TRI(i,:)];
    end
end
TRI = TRINOW;


t1 = find(abs(Nodes(:,1) + XX) <1e-10);
t2 = find(abs(Nodes(:,1)-(XX+1)) <1e-10);

DirichletNodes = [t1; t2];
%===
% Figures
%===



plot(OuterBoundary(:,1), OuterBoundary(:,2),'-*b')
hold on
plot(InnerBoundary(:,1), InnerBoundary(:,2),'-or')

plot(Nodes(:,1), Nodes(:,2),'g.')

axis equal

figure
trisurf(TRI, Nodes(:,1), Nodes(:,2), Nodes(:,1)*0)
axis equal
view([0 0 1])