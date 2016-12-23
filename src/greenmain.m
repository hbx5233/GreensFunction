%--------------------------------------------------------------------------

% Copyright University of Chinese Academy of Sciences,2015.
% All rights reserved.
% greenmain.m : $Date:2015/04/29 18:48:30$

%	Credits:UCAS Lin Yi, April 29th,2015
%	Reference:
%  Lane R. Johnson,
%  1974,Green's Function for Lamb's Problem,Geophys,J.R.Soc.37,99-131.

clear;
tic;            %start timing

%----------------------------parameter setting-----------------------------
%----------------------------unchange parameter----------------------------
x3=0;           %record x3 direction coordinate
x10=0;          %source x1 direction coordinate
x20=0;          %source x2 direction coordinate
t0=0;           %initial time

rho=3.3;        %density
alpha=8.0;      %P-wave velocity
beta=4.62;      %S-wave velocity
mu=rho*beta^2;  % shear modulus
F=[1;1;1];      %body force value,set as a unit force,not used in
                %this program

dt=(4-1)/300;   %time sample interval
t=1:dt:4;       % time value according to the fig 2-4 

%---------------------------input parameter--------------------------------
%-------------------calculate g(2,0,0,t;0,0,10,0)--------------------------
x1=2;           %record x1 direction coordinate
x2=0;           %record x2 direction coordinate
x30=10;         %source x3 direction coordinate

%--------parameter calculated by the previous parameter--------------------
r=sqrt(x1^2+x2^2+(x3-x30)^2);   %distance from source to record point
theta=acos((x30-x3)/r);         %included angle between vetor r and x3 axis
phi=atan(x2/x1);                %angle made by projection of record point 
                                %and x1 axis

g=calgreen(r,alpha,beta,theta,phi,mu,F,dt);       %calculate green function
interval=[11 5 7 9 1]*1e-4;                       %set plot interval
%--------------loop to update the g_(ij) by add plot interval--------------
for i=1:5
    g(i,:)=interval(i)+g(i,:);
end
%-----------------------plot figure 2 in the paper-------------------------
figure;
plot(t,g(1,:));
text(0.6,interval(1),'g_{11}^H');
hold on;
plot(t,g(2,:));
text(0.6,interval(2),'g_{13}^H');
plot(t,g(3,:));
text(0.6,interval(3),'g_{22}^H');
plot(t,g(4,:))
text(0.6,interval(4),'g_{31}^H');
plot(t,g(5,:))
text(0.6,interval(5),'g_{33}^H');
axis([0.5 4 0 1.4e-3]);
set(gca,'XTick',0.5:0.5:4, ...
        'YTick',0:1e-4:1.4e-3); 
set(gca,'xticklabel',{'','1','','2','','3','','4'})
set(gca,'yticklabel','');
xlabel('t,sec');  
%title('The component of G^H(2,0,0,t;0,0,10,0)'); 
title('Johnson(1974)'); 
print('Johnson_fig2','-dpdf');

%-------------------------calculate g(10,0,0,t;0,0,2,0)--------------------
x1=10;
x2=0;
x30=2;

%--------------parameter calculated by the previous parameter--------------
r=sqrt(x1^2+x2^2+(x3-x30)^2);
theta=acos((x30-x3)/r);
phi=atan(x2/x1);


g=calgreen(r,alpha,beta,theta,phi,mu,F,dt);
interval=[11 6 8 10 2]*1e-4;
for i=1:5
    g(i,:)=interval(i)+g(i,:);
end
%--------------------------plot figure 3 in the paper----------------------
figure;
plot(t,g(1,:));
text(0.6,interval(1),'g_{11}^H');
hold on;
plot(t,g(2,:));
text(0.6,interval(2),'g_{13}^H');
plot(t,g(3,:));
text(0.6,interval(3),'g_{22}^H');
plot(t,g(4,:))
text(0.6,interval(4),'g_{31}^H');
plot(t,g(5,:))
text(0.6,interval(5),'g_{33}^H');
axis([0.5 4 0 1.4e-3]);
set(gca,'XTick',0.5:0.5:4, ...
        'YTick',0:1e-4:1.4e-3); 
set(gca,'xticklabel',{'','1','','2','','3','','4'})
set(gca,'yticklabel','');
xlabel('t,sec');  
%title('The component of G^H(10,0,0,t;0,0,2,0)'); 
title('Johnson(1974)'); 
print('Johnson_fig3','-dpdf');

%----------------------calculate g(10,0,0,t;0,0,0.2,0)---------------------
x1=10;
x2=0;
x30=0.2;

%-----------parameter calculated by the previous parameter-----------------
r=sqrt(x1^2+x2^2+(x3-x30)^2);
theta=acos((x30-x3)/r);
phi=atan(x2/x1);


g=calgreen(r,alpha,beta,theta,phi,mu,F,dt);
interval=[16 12 13 15 6]*1e-4;
for i=1:5
    g(i,:)=interval(i)+g(i,:);
end
%-----------------------plot figure 4 in the paper-------------------------
figure;
plot(t,g(1,:));
text(0.6,interval(1),'g_{11}^H');
hold on;
plot(t,g(2,:));
text(0.6,interval(2),'g_{13}^H');
plot(t,g(3,:));
text(0.6,interval(3),'g_{22}^H');
plot(t,g(4,:))
text(0.6,interval(4),'g_{31}^H');
plot(t,g(5,:))
text(0.6,interval(5),'g_{33}^H');
axis([0.5 4 0 1.9e-3]);
set(gca,'XTick',0.5:0.5:4, ...
        'YTick',0:1e-4:1.9e-3); 
set(gca,'xticklabel',{'','1','','2','','3','','4'})
set(gca,'yticklabel','');
xlabel('t,sec');  
%title('The component of G^H(10,0,0,t;0,0,0.2,0)'); 
title('Johnson(1974)'); 
print('Johnson_fig4','-dpdf');

tt1=toc;                                           %end of timing
display(strcat('total time£º',num2str(tt1),'sec'));%show total running time
                                                   %command window

