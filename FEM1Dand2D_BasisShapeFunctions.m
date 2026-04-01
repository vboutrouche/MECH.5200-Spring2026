clear all 
close all

xi = -1:0.05:1;
eta = -1:0.05:1;

%=========================================================================
% 1D Basis Functions -- Calculating
%=========================================================================

% Let's calculate and plot the two linear basis functions
N1 = 1/2*(1-xi)
N2 = 1/2*(1+xi)

NQ1 = xi/2.*(xi-1);
NQ2 = (1-xi.^2);
NQ3 = xi/2.*(xi+1);

%=========================================================================
% 1D Basis Functions -- Plotting 
%=========================================================================
figure(1)
subplot(1,2,1)
plot(xi, N1,'r-','linewidth',2)
hold on
plot(xi, N2,'b-','linewidth',2)
plot(xi, N1*2+N2*3,'m-','linewidth',3)
plot(1,5,'w.')
plot(1,-1,'w.')
plot([-1 1],[0 0],'-o', 'color',[.25 .25 .25],'linewidth',3)
grid on
xlabel('\xi')
title('Linear Basis Functions')
a = legend('N1','N2','\alpha \cdot N_1 + \beta \cdot N_2 ')
set(a,'Location','northwest')
set(gca,'fontsize',16)

subplot(1,2,2)
plot(xi, NQ1,'r-','linewidth',2)
hold on
plot(xi, NQ2,'b-','linewidth',2)
plot(xi, NQ3,'m-','linewidth',2)
plot(xi, NQ1*2+NQ2*2.75+NQ3*3,'g-','linewidth',3)
plot(1,5,'w.')
plot(1,-1,'w.')
plot([-1 1],[0 0],'o-', 'color',[.25 .25 .25],'linewidth',3)
grid on
xlabel('\xi')

title('Quadratic Basis Functions')
a= legend('N1','N2', 'N_3','\alpha \cdot N_1 + \beta \cdot N_2 + \gamma \cdot N_3')
set(a,'Location','northwest')
set(gca,'fontsize',16)

%=========================================================================
% 2D Linear Basis Functions -- Calculating
%=========================================================================

for(i = 1:length(xi))
    for(j = 1:length(eta))
        N2_1(i,j) = 1/4*(1-xi(i)).*(1-eta(j));
        N2_2(i,j) = 1/4*(1+xi(i)).*(1-eta(j));
        N2_3(i,j) = 1/4*(1+xi(i)).*(1+eta(j));
        N2_4(i,j) = 1/4*(1-xi(i)).*(1+eta(j));
        XI(i,j) = xi(i);
        ETA(i,j) = eta(j);
    end
end


%=========================================================================
% 2D Linear Basis Functions -- Plotting
%=========================================================================
figure
subplot(2,2,1)

%--------------------------
% Basis 1
%--------------------------
%plot3(XI, ETA, N2_1,'*b')
a = surf(XI, ETA, N2_1)
set(a,'facecolor','b','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_1')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light

%--------------------------
% Basis 2
%--------------------------
subplot(2,2,2)
% plot3(XI, ETA, N2_2,'*r')
a = surf(XI, ETA, N2_2)
set(a,'facecolor','r','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_2')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light

%--------------------------
% Basis 3
%--------------------------
subplot(2,2,3)
%plot3(XI, ETA, N2_3,'*m')
a = surf(XI, ETA, N2_3)
set(a,'facecolor','m','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_3')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light

%--------------------------
% Basis 4
%--------------------------
subplot(2,2,4)
% plot3(XI, ETA, N2_4,'*g')
a = surf(XI, ETA, N2_4)
set(a,'facecolor','g','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_4')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light

%--------------------------
% Weighted Sum of the basis functions
%--------------------------
figure
% plot3(XI, ETA, N2_1*2+N2_2*3+N2_3*1+N2_4*1.5,'*g')
a = surf(XI, ETA, N2_1*2+N2_2*3+N2_3*1+N2_4*1.5)
set(a,'facecolor','g','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('\Sigma w_i*N_i = 2\cdot N_1 + 3\cdot N_2 + N_3 + 1.5\cdot N_4')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light



%=========================================================================
% 2D Quadratic Basis Functions -- Calcualting (9-point)
%=========================================================================
for(i = 1:length(xi))
    for(j = 1:length(eta))
        N2q_1(i,j) = 1/4*xi(i).*eta(j).*(xi(i)-1).*(eta(j)-1);
        N2q_2(i,j) = eta(j)/2.*(1-xi(i).^2).*(eta(j)-1);
        N2q_3(i,j) = (1-xi(i).^2).*(1-eta(j).^2);
        
        XI(i,j) = xi(i);
        ETA(i,j) = eta(j);
    end
end


%=========================================================================
% 2D Quadratic Basis Functions -- Plotting (9-point)
%=========================================================================



figure

%--------------------------
% Basis 1
%--------------------------
subplot(2,2,1)
%plot3(XI, ETA, N2q_1,'*b')
a = surf(XI, ETA, N2q_1)
set(a,'facecolor','b','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_1')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light


%--------------------------
% Basis 2
%--------------------------
subplot(2,2,2)
%plot3(XI, ETA, N2q_2,'*r')
a = surf(XI, ETA, N2q_2)
set(a,'facecolor','r','facealpha',0.75)

hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_2')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light


%--------------------------
% Basis 3
%--------------------------
subplot(2,2,3)
% plot3(XI, ETA, N2q_3,'*m')
a = surf(XI, ETA, N2q_3)
set(a,'facecolor','m','facealpha',0.75)

hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_3')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light


%--------------------------
% Basis weighted combination
%--------------------------
subplot(2,2,4)
% plot3(XI, ETA, N2q_1*1.5 + N2q_2*1 + N2q_3*2,'*g')
a = surf(XI, ETA, N2q_1*1.5 + N2q_2*1 + N2q_3*2)
set(a,'facecolor','g','facealpha',0.75)

hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('1.5\cdot N_1 + 1.0\cdot N_2 + 2.0\cdot N_3')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light
%================================








%=========================================================================
% 2D Linear Basis Functions -- Plotting
%=========================================================================
figure


%--------------------------
% Basis 1
%--------------------------
%plot3(XI, ETA, N2_1,'*b')
a = surf(1+XI, 1+ETA, N2_1)
set(a,'facecolor','g','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1]+1,[-1 1 1 -1 -1]+1,[0 0 0 0 0],'k');
%xset(a,'facealpha',0.65)
xlabel('\xi')
ylabel('\eta')
title('Basis N_1')
set(gca,'fontsize',16)
view([1. 1.25 .75])


%--------------------------
% Basis 2
%--------------------------
hold on
% plot3(XI, ETA, N2_2,'*r')
a = surf(1+XI-2, 1+ETA, N2_2)
set(a,'facecolor','g','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1]+1-2,1+[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.65)
xlabel('\xi')
ylabel('\eta')
title('Basis N_2')
set(gca,'fontsize',16)
view([1. 1.25 .75])


%--------------------------
% Basis 3
%--------------------------
%plot3(XI, ETA, N2_3,'*m')
a = surf(1+XI-2, 1+ETA-2, N2_3)
set(a,'facecolor','g','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1]+1-2,1-2+[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.65)
xlabel('\xi')
ylabel('\eta')
title('Basis N_3')
set(gca,'fontsize',16)
view([1. 1.25 .75])


%--------------------------
% Basis 4
%--------------------------

% plot3(XI, ETA, N2_4,'*g')
a = surf(1+XI, ETA+1-2, N2_4)
set(a,'facecolor','g','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1]+1,1-2+[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.65)
xlabel('\xi')
ylabel('\eta')
title('Basis N_4')
set(gca,'fontsize',16)
view([1. 1.25 .75])



%--------------------------
% Basis 1
%--------------------------
%plot3(XI, ETA, N2_1,'*b')
a = surf(1+XI+1, 1+ETA+1, N2_1)
set(a,'facecolor','r','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1]+2,2+[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_1')
set(gca,'fontsize',16)
view([1. 1.25 .75])


%--------------------------
% Basis 2
%--------------------------
hold on
% plot3(XI, ETA, N2_2,'*r')
a = surf(1+XI-2+1, 1+ETA+1, N2_2)
set(a,'facecolor','r','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],2+[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_2')
set(gca,'fontsize',16)
view([1. 1.25 .75])


%--------------------------
% Basis 3
%--------------------------
%plot3(XI, ETA, N2_3,'*m')
a = surf(1+XI-2+1, 1+ETA-2+1, N2_3)
set(a,'facecolor','r','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1],[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_3')
set(gca,'fontsize',16)
view([1. 1.25 .75])


%--------------------------
% Basis 4
%--------------------------

% plot3(XI, ETA, N2_4,'*g')
a = surf(1+XI+1, ETA+1-2+1, N2_4)
set(a,'facecolor','r','facealpha',0.75)
hold on
a = patch([1 1 -1 -1 1]+2,[-1 1 1 -1 -1],[0 0 0 0 0],'k');
%set(a,'facealpha',0.5)
xlabel('\xi')
ylabel('\eta')
title('Basis N_4')
set(gca,'fontsize',16)
view([1. 1.25 .75])
light
