function F=F_angle(phi,theta)
global Size_X Size_Y E X Y k
% F=0;
% if (length(theta)>1)
%     aaa=1;
% end
% for i=1:Size_X
%     for j=1:Size_Y
%         F=F+E(i,j)*exp(1i*k*(X(i,j)*cos(phi)+Y(i,j)*sin(phi))*sin(theta))*sin(theta); %% sin(theta) вне квадрата модуля!!
%     end
% end
% F=(abs(F)).^2;
for i=1:length(phi)
    F(1,i)=( abs( sum( sum( E.*exp(1i.*k.*( X.*cos(phi(i)) + Y.*sin(phi(i))).*sin(theta)) ) ) ) )^2.*sin(theta);
end