%This is a function of the integrand of the equation(26)
function integf=integrand(v,t,r,alpha,beta,theta,phi,F,MN,index_int)
%self document
%v    :integral variable
%t    :time
%r    :distance from source to record point
%alpha:P-wave velocity
%beta :S-wave velocity
%theta:included angle between vetor r and x3 axis
%phi  :angle made by projection of record point and x1 axis
%mu   :shear modulus
%F    :body force value,set as a unit force,not used in this program
%MN   :{'M11','M13','M22','M31','M33',...
       %'N11','N13','N22','N31','N33'} 
%index_int: %index of integral corresponding to the four 
            %integral in equation(26) respectly.

switch(index_int);%transform the variable of integral in the first,
                  %second and forth integral in equation(26)
  case '1';
    p=sqrt((t/r)^2-alpha^-2)-v.^2;
  case '2';
    p=sqrt((t/r)^2-beta^-2)-v.^2;
  case '3';
    p=v;          %do not need
  case '4';
    p=sqrt((t/r)^2-beta^-2)+v.^2;
end

%-------calculate the different value of q in different integral-----------
switch(index_int)
  case '1';
    q=-t/r*sin(theta)+(sqrt((t/r)^2-alpha^-2-p.^2)*cos(theta))*i;
  case '2';
    q=-t/r*sin(theta)+(sqrt((t/r)^2-beta^-2-p.^2)*cos(theta))*i;
  case '3';
    q=-t/r*sin(theta)+(sqrt(-(t/r)^2+beta^-2+p.^2)*cos(theta));
  case '4';
    q=-t/r*sin(theta)+(sqrt(-(t/r)^2+beta^-2+p.^2)*cos(theta));
end

eta_alpha=sqrt(alpha^-2+p.^2-q.^2);
eta_beta=sqrt(beta^-2+p.^2-q.^2);
gamma=eta_beta.^2+p.^2-q.^2;
sigma=gamma.^2+4*eta_alpha.*eta_beta.*(q.^2-p.^2);

switch(MN);
  case 'M11';
    M=2*eta_beta.*((q.^2+p.^2)*(cos(phi))^2-p.^2);
  case 'M12';
    M=2*eta_beta.*(q.^2+p.^2)*sin(phi)*cos(phi);
  case 'M21';
    M=2*eta_beta.*(q.^2+p.^2)*sin(phi)*cos(phi);
  case 'M13';
    M=2*q.*eta_alpha.*eta_beta*cos(phi);
  case 'M22';
    M=2*eta_beta.*((q.^2+p.^2)*(sin(phi))^2-p.^2);
  case 'M23';
    M=2*q.*eta_alpha.*eta_beta.*sin(phi);
  case 'M31';
    M=q.*gamma*cos(phi);
  case 'M32';
    M=q.*gamma*sin(phi);
  case 'M33';
    M=eta_alpha.*gamma;
  case 'N11';
    N=(eta_beta.^2.*gamma-(gamma-4*eta_alpha.*eta_beta).*((q.^2+p.^2).*...
        (sin(phi))^2-p.^2))./eta_beta;
  case 'N12';
    N=(q.^2+p.^2).*(gamma-4*eta_alpha.*eta_beta)*sin(phi)*cos(phi)./...
        eta_beta;
  case 'N21';
    N=(q.^2+p.^2).*(gamma-4*eta_alpha.*eta_beta)*sin(phi)*cos(phi)./...
        eta_beta;
  case 'N13';
    N=-q.*gamma*cos(phi);
  case 'N22';
    N=(eta_beta.^2.*gamma-(gamma-4*eta_alpha.*eta_beta).*((q.^2+p.^2).*...
        (cos(phi))^2-p.^2))./eta_beta;
  case 'N23';
    N=-q.*gamma*sin(phi);
  case 'N31';
    N=-2*q.*eta_alpha.*eta_beta*cos(phi);
  case 'N32';
    N=-2*q.*eta_alpha.*eta_beta*sin(phi);
  case 'N33';
    N=2*eta_alpha.*(q.^2-p.^2);
end

%calculate different integral 
switch(index_int);
  case '1';
    SRR=real(M.*eta_alpha.*sigma.^-1.*(2*sqrt((t/r)^2-alpha^-2)-v.^2).^-0.5);
    integf=SRR*heaviside(t-r/alpha);
  case '2';
    SRR=real(N.*eta_beta.*sigma.^-1.*(2*sqrt((t/r)^2-beta^-2)-v.^2).^-0.5);
    integf=SRR*heaviside(t-r/beta);
  case '3';
    t2=r/alpha*sin(theta)+r*sqrt(beta^-2-alpha^-2)*cos(theta);
    SRR=imag(N.*eta_beta.*sigma.^-1.*(beta^-2+p.^2-(t/r)^2).^-0.5);
    integf=SRR*heaviside(sin(theta)-beta/alpha)*(heaviside(t-t2)-...
       heaviside(t-r/beta));
  case '4';
    SRR=imag(N.*eta_beta.*sigma.^-1.*(2*sqrt((t/r)^2-beta^-2)+v.^2).^-0.5);
    integf=SRR*heaviside(sin(theta)-beta/alpha)*heaviside(t-r/beta);
end

 
 end
