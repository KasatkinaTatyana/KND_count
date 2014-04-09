close all
clear all
clc
global lambda N1 N2 dx dy
% lambda=3.37;
f = 34e9;       %[��]
c = 3e11;       %[��/�]
lambda = c/f;   %[��]
N1=50;
N2=50;
dx=0.7*lambda;
dy=0.6*lambda;
%----------------------������� ��������------------------------------------
% tic;
% I=dblquad(@Int2Sum,0,2*pi,0,pi/2);
% toc;
% disp('Double integral: ');
% disp(I);
%---------------------------�����------------------------------------------
tic;
S=0;
S_max=0;
%--------------------------------------------------------------------------
A1=1;
A2=1;
%--------------------------------------------------------------------------
x_c=dx*(N1-1)/2;
y_c=dy*(N2-1)/2;
r_max=(x_c^2+y_c^2)^0.5;
count=0;
%--------------------------------------------------------------------------
for x1=0:dx:(N1-1)*dx
    for y1=0:dy:(N2-1)*dy
        for x2=0:dx:(N1-1)*dx
            for y2=0:dy:(N2-1)*dy
                x=abs(x1-x2);
                y=abs(y1-y2);
                %----------------------------------------------------------
%                 r1=((x1-x_c)^2+(y1-y_c)^2)^0.5;
%                 r2=((x2-x_c)^2+(y2-y_c)^2)^0.5;
%                 
%                 A1=cos(pi/2-pi/2*r1/r_max);
%                 A2=cos(pi/2-pi/2*r2/r_max);
                %----------------------------------------------------------
                qxy=(x^2+y^2)^0.5;
                               
                if ((abs(x)<dx/10)&&(abs(y)<dy/10))
                    S=S+A1*A2;
                    count=count+1;
                else
                    S=S+A1*A2*sin(qxy)/qxy;
                end
            end
        end
        S_max=S_max+A1^2;
%         S=S+A1^2;
    end
end                                                                                                                                          
KND=2*S_max^2/S
% KND=2*pi*S_max^2/S   ??? 
KND_log=10*log10(KND)
toc;
% disp('Summa:');
% disp(2*pi*S);
