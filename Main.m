close all
clear all
clc
global lambda N1 N2 dx dy
% lambda=3.37;
f = 34e9;       %[Гц]
c = 3e11;       %[мм/с]
lambda = c/f;   %[мм]
% k = 2*pi/lambda;
k = 1;

dx=0.7*lambda;
dy=0.6*lambda;

Size_X=2;
Size_Y=2;

A=dx*Size_X;   
B=dy*Size_Y;    

N1 = Size_X+1;
N2 = Size_Y+1;

tic;

for count1=1:N1
    for count2=1:N2
        X(count1,count2) = (-dx*((N1-1)/2) + dx*(count1-1));
        Y(count1,count2) = (dy*((N2-1)/2) - dy*(count2-1));
    end
end

% figure
% plot(X,Y,'*r')
%% Определение фазы


%theta_otkl = 40;
theta_otkl = 0;
phi_otkl = 0;

% figure
% plot3(X,Y,E,'*r')

Phase = - k.* (X.*cos(pi/180*phi_otkl) + Y.*sin(pi/180*phi_otkl)).*sin(pi/180*theta_otkl);
ex1 = exp(1i*Phase);

dtheta = 0.1;
theta = 0:dtheta:90;
dphi = 1;
phi = 0:dphi:360;

% Fizl = cos(theta*pi/180);
Fizl=1;
%----------------------Двойной интеграл------------------------------------
% tic;
% I=dblquad(@Int2Sum,0,2*pi,0,pi/2);
% toc;
% disp('Double integral: ');
% disp(I);
%--------------------------------------------------------------------------
x_c=dx*(N1-1)/2;
y_c=dy*(N2-1)/2;

%Определение амплитуды
% E1 = zeros(N1,N2);
% E2 = zeros(N1,N2);
% for count1=1:N2
%     E1(:,count1) = 0.3 + 0.7*cos(X(:,count1)/x_c*pi/2);
% end
% for count1=1:N1
%     E2(count1,:) = 0.3 + 0.7*cos(Y(count1,:)/y_c*pi/2);
% end
% 
% E = E1.*E2;
E = ones(N1,N2);

%%КНД
%---------------------------Сумма------------------------------------------
S=0;
S_max=0;
%--------------------------------------------------------------------------
for c1=1:N1
    for c2=1:N2
        for c3=1:N1
            for c4=1:N2
                x=abs(X(c1,c2)-X(c3,c4));
                y=abs(Y(c1,c2)-Y(c3,c4));
                %----------------------------------------------------------
                qxy=k*(x^2+y^2)^0.5;
                   
                S=S+E(c1,c2)*E(c3,c4)*sinc(qxy/pi);
            end
        end
        S_max=S_max+E(c1,c2);       
    end
end                                                                                                                                          
G=2*S_max^2/S

G_log=10*log10(G)
%% Вычисление ДН и КНД(theta,phi)

for count2 = 1:length(phi)
    FS = zeros(1,length(theta));
    for count1 = 1:length(theta)
        Exp = exp(1i.*k.*( X.*cos(phi(count2).*pi./180) + Y.*sin(phi(count2).*pi./180)).*sin(theta(count1).*pi./180) );
        FS(count1) = sum(sum(E.*ex1.*Exp));
    end
    F{count2} = FS.*Fizl;
    F1{count2} = (abs(F{count2})./S_max);     % Нормированная ДН
    G_angle{count2} = 2*((abs(F1{count2})).^2)./S;
end

toc;
%% ------------------------Проверка-----------------------------------------
% F_max=F1{1}(1);
% for count1=1:length(phi)
%     for count2=1:length(theta)
%         G_angle_2{count1}(count2)=2*(abs(F1{count1}(count2)))^2/S;
%         differ = G_angle_2{count1}(count2) - G_angle{count1}(count2);
%         if (differ ~= 0)
%             disp(differ);
%         end
%         
%         if F_max<F1{count1}(count2)
%             F_max=F1{count1}(count2);
%         end
%     end
% end
% 
% F_max

% disp('Summa:');
% disp(2*pi*S);

