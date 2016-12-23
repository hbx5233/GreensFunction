%This function is to calculate the green's function for a source which is a
%unit step function in time

function gg=calgreen(r,alpha,beta,theta,phi,mu,F,dt)

%self document
%r    :distance from source to record point
%alpha:P-wave velocity
%beta :S-wave velocity
%theta:included angle between vetor r and x3 axis
%phi  :angle made by projection of record point and x1 axis
%mu   :shear modulus
%F    :body force value,set as a unit force,not used in this program
%dt   :time sample interval                

t=1:dt:4;                     %separate time by dt
%------------------------zero out integral array---------------------------
%INTEGF1,INTEGF2,INTEGF3,INTEGF4 corresponds to the four integral in 
%equation(26) respectly.

INTEGF1=zeros(length(t),5);
INTEGF2=zeros(length(t),5);
INTEGF3=zeros(length(t),5);
INTEGF4=zeros(length(t),5);

MN={'M11','M13','M22','M31','M33',...
    'N11','N13','N22','N31','N33'}; %same meaning as in the paper
                                                 
index_int={'1','2','3','4'};  %index of integral corresponding to the four 
                              %integral in equation(26) respectly.
                              
%calculate the first integral in equation(26), P-wave
  for m=1:1:length(t)
      limup=((t(m)/r)^2-alpha^-2)^0.25; %transformation of the
                                        %limits of integration
                                        %according to equation(134)
%use quadl method to calculate the integral
%integrand() is the integrand function
%quadl(fun(),a,b),fun() is the integrand function
%a,b stand for a lower or upper limit of integration respectly


      for i=1:5
          INTEGF1(m,i)=quadl(@(v)integrand(v,t(m),r,alpha,beta,theta,...
          phi,F,MN{i},index_int{1}),0,limup);
      end 
  end
  


%calculate the second integral in equation(26), 
for m=1:1:length(t)
    limup2=((t(m)/r)^2-beta^-2)^0.25;%transform like the first integral
	for i=1:5
        INTEGF2(m,i)=quadl(@(v)integrand(v,t(m),r,alpha,beta,theta,phi,F,...
        MN{i+5},index_int{2}),0,limup2);
    end
end

%calculate the third integral in equation(26), 
for m=1:1:length(t)
	p2(m)=sqrt(((t(m)/r-sqrt(beta^-2-alpha^-2)*cos(theta))/sin(theta))^2 ...
          -alpha^-2);
    for i=1:5
      INTEGF3(m,i)=quadl(@(v)integrand(v,t(m),r,alpha,beta,theta,phi,...
                   F,MN{i+5},index_int{3}),0,p2(m));
    end
    
end

%calculate the forth integral in equation(26), 
for m=1:1:length(t)
  p2=sqrt(((t(m)/r-sqrt(beta^-2-alpha^-2)*cos(theta))/sin(theta))^2 ...
      -alpha^-2);
  p1=sqrt((t(m)/r)^2-beta^-2);
  limup4=sqrt(p2-p1);
 for i=1:5
    INTEGF4(m,i)=quadl(@(v)integrand(v,t(m),r,alpha,beta,theta,phi,F,...
        MN{i+5},index_int{4}),0,limup4);
 end
end

S=1/(pi^2*mu*r);       %constant index in equation(26)
% add up the four integrals to get the green's function
for m=1:length(t)
    gg(:,m)=S*(2*INTEGF1(m,:)+2*INTEGF2(m,:)-INTEGF3(m,:)-2*INTEGF4(m,:));
    %gg is the green's function;
    %the factor 2 in this equation  because we perform a
    %transformation of the variable of integration
end
end

