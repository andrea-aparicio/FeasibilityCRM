%% colors 

colors = ["#BC4447", "#6D9BE2", "#2b9f4d", "#542d89", "#b5b443", "#9950d7", "#6f8f2b", "#c46029", "#cb4eb7", "#5ec1cd"];
lw=1;
%% 1
At=readmatrix('Asimst.csv','NumHeaderLines',1);
A=At(:,2:11);

%N = .9.*rand(10,1)+.1;
N = [0.8566,
    0.3289,
    0.8329,
    0.3192,
    0.9363,
    0.4150,
    0.2769,
    0.3260,
    0.6544,
    0.5260]
r = -A*N;



%x0=unifrnd(.1,.5,[10,1]);
x0 = [    0.2407,
    0.4323,
    0.3341,
    0.3199,
    0.4669,
    0.2143,
    0.4029,
    0.4015,
    0.2522,
    0.3271]


params.A=A;
params.r=r;
f=@(t,x)MACRFunc(x,params);

t0=0;
tf=1000000;
dt = 1;

[x,t]=RK4(f,x0,t0,tf,dt);
figure()
loglog(t,x)
%% 2

A2 = A;
nrand = abs(normrnd(0,.3,5,1));
for i=1:5
    A2(i,i)=-nrand(i);
end

N2 = N;% .9.*rand(10,1)+.1
r2 = -A2*N2
params2.A=A2;
params2.r=r2;
f2=@(t,x2)MACRFunc(x2,params2);


[x2,t]=RK4(f2,x0,t0,tf,dt);
figure()
loglog(t,x2)
%% PLOTS 1
xp=x;
lw=1;
subplot(2,2,3)

loglog(t,xp(1,:),'Color',colors(1), 'LineWidth',lw)
hold
for i=2:5
    loglog(t,xp(i,:),'Color',colors(i), 'LineWidth',lw)
   %semilogx(t,x(i,:),'Color',colors(i), 'LineWidth',lw)
    
end
xlim([0,10E5])
ylim([10E-3,10E0])
xticks([0,1, 10E1, 10E2, 10E3, 10E4])
yticks([10E-3,10E-2,10E-1,10E0, 10E1])
xlabel('time')
ylabel({'resource concentration'; 'mass/vol.'})
set(gca,'FontSize',14)
hold


subplot(2,2,1)
loglog(t,xp(6,:),'Color',colors(6), 'LineWidth',lw)
hold
for i=7:10
    loglog(t,xp(i,:),'Color',colors(i), 'LineWidth',lw)
    
end
ylim([5E-2,2E0])
xlim([0,10E4])
ylim([10E-2,2])
xticks([0,1, 10E1, 10E2, 10E3, 10E4])
yticks([10E-3,10E-2,10E-1,10E0, 10E1])
hold
ylim([-.01,1])
x8_2=x;
xlabel('time')
ylabel({'species density'; 'indiv./vol.'})
set(gca,'FontSize',14)
min(min(x))
% writematrix(x,'x40.3.csv')

%% Plots2

xp=x2;
lw=1;
subplot(2,2,4)
loglog(t,xp(1,:),'Color',colors(1), 'LineWidth',lw)
hold
for i=2:5
    loglog(t,xp(i,:),'Color',colors(i), 'LineWidth',lw)
   
    
end
xlim([0,10E4])
ylim([40E-3,5E0])
xticks([0,1, 10E1, 10E2, 10E3, 10E4])
yticks([10E-3,10E-2,10E-1,10E0, 10E1])
xlabel('time')
ylabel({'resource concentration'; 'mass/vol.'})
set(gca,'FontSize',14)
hold


subplot(2,2,2)
loglog(t,xp(6,:),'Color',colors(6), 'LineWidth',lw)
hold
for i=7:10
    loglog(t,xp(i,:),'Color',colors(i), 'LineWidth',lw)
    
end
ylim([5E-2,2E0])
xlim([0,10E4])
ylim([10E-2,2])

xticks([0,1, 10E1, 10E2, 10E3, 10E4])
yticks([10E-3,10E-2,10E-1,10E0, 10E1])
hold
%ylim([-.01,1])
x8_2=x;
xlabel('time')
ylabel({'species density'; 'indiv./vol.'})
set(gca,'FontSize',14)
min(min(x))
% writematrix(x,'x40.3.csv')

