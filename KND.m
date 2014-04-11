clc
clear all
close all
global Size_X Size_Y E X Y k

f = 34e9;       %[Гц]
c = 3e11;       %[мм/с]
lambda = c/f;   %[мм]

% k = 2*pi/lambda;

k = 1;
 
% dx = 6;         %[мм]
% dy = 6;         %[мм]
% 
% A = 60;
% B = 60;

dx=lambda*0.7;
dy=lambda*0.6;

Size_X=2;
Size_Y=2;
A=dx*Size_X;    % Size_X = N1 - 1
B=dy*Size_Y;    % Size_Y = N2 - 1

m = floor(A/dx);
n = floor(B/dy);

tic
for count1=1:m+1
    for count2=1:n+1
        X(count1,count2) = (-dx*(m/2) + dx*(count1-1));
        Y(count1,count2) = (dy*(n/2) - dy*(count2-1));
    end
end

m = m+1;
n = n+1;

% figure
% plot(X,Y,'*r')
%% Определение фазы и амплитуды


%theta_otkl = 40;
theta_otkl = 0;
phi_otkl = 0;



% E1 = zeros(m,n);
% E2 = zeros(m,n);
% for count1=1:n
%     E1(:,count1) = 0.3 + 0.7*cos(X(:,count1)/abs(min(min(X)))*pi/2);
% end
% for count1=1:m
%     E2(count1,:) = 0.3 + 0.7*cos(Y(count1,:)/abs(min(min(Y)))*pi/2);
% end
% 
% E = E1.*E2;
%% 
% figure
% plot3(X,Y,E,'*r')
%%

E = ones(m,n);

Phase = - k.* (X.*cos(pi/180*phi_otkl) + Y.*sin(pi/180*phi_otkl)).*sin(pi/180*theta_otkl);
ex1 = exp(1i*Phase);

dtheta = 0.1;
theta = 0:dtheta:90;
dphi = 1;
phi = 0:dphi:360;

% Fizl = cos(theta*pi/180);
Fizl=1;

%% Определение максимума ДН

FN = zeros(1,length(theta));
for count1 = 1:length(theta)
    Exp = exp(1i.*k.*( X.*cos(phi_otkl.*pi./180) + Y.*sin(phi_otkl.*pi./180)).*sin(theta(count1).*pi./180) );
%     FN(count1) = sum(sum(E.*ex1.*Exp));
    FN(count1) = sum(sum(E));
end

Norm = max(abs(FN.*Fizl));


%% Вычисление ДН

for count2 = 1:length(phi)
    FS = zeros(1,length(theta));
    for count1 = 1:length(theta)
        Exp = exp(1i.*k.*( X.*cos(phi(count2).*pi./180) + Y.*sin(phi(count2).*pi./180)).*sin(theta(count1).*pi./180) );
        FS(count1) = sum(sum(E.*ex1.*Exp));
    end
    F{count2} = FS.*Fizl;
%     F1{count2} = (abs(F{count2})./Norm);
    F1{count2} = (abs(F{count2})./1);
%     F3{count2} = 20.*log10(F1{count2});
end

%% Расчет КНД


% D = 0;
% for count1=1:length(phi)
%     D = D + sum(F1{count1}.^2 .*sin(theta*pi/180))*dtheta*dphi*pi*pi/180/180;
% end
D=dblquad(@F_angle,0,2*pi,0,pi/2,0.001);

D1 = 4*pi*Norm^2/D;

disp(D1);
D_dB = 10*log10(D1)
%% Объемный КНД
for count2 = 1:length(phi)
    F1{count2} = (abs(F{count2})./Norm);
    G_angle{count2} = 4*pi*((abs(F1{count2})).^2)./D;
end
toc;
%% Адекватность результата

D2 = 4*pi*A*B/lambda/lambda;
D2_dB = 10*log10(D2)


