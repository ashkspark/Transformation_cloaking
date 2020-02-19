clear all; clc;
mu = 1.0;
lam = 1.0;
rho = 1.0;
omega = 30.00;
am = 20.0;
Rt = 2.0;
Rs = 0.7*Rt;
epsi = Rs/500; %1/100 of the finite hole
m = 400;
delta = 0.025*Rt;
Ri = Rs/5;
Ro = 2*Ri;
del = Ro-epsi;
k = (rho*omega^2/mu)^(1/2);
s = 191;
g=(s+1)/2;
[r,t] = meshgrid(epsi:(Rt-epsi)/m:Rt,-pi:2*pi/m:pi);
z = ones(length(r),length(t));
x = r.*cos(t);
y = r.*sin(t);
z((x-Rs).^2+y.^2<delta.^2) = NaN;
x(isnan(z)) = NaN;
y(isnan(z)) = NaN;
for n = 1:1:s
Uii = am./4./1i./mu.*besselh(n-g,k.*Rs).*besselj(n-g,k.*r).*exp(1i.*(n-g).*t);
Uii(r > Rs) = 0;
Uio = am./4./1i./mu.*besselh(n-g,k.*r).*besselj(n-g,k.*Rs).*exp(1i.*(n-g).*t);
Uio(r <= Rs) = 0;
Ui = Uii+Uio;

Us = (-1i).^(n-g).*besselh(n-g,k.*r).*exp(1i.*(n-g).*t)...
    .*((-1).^(n-g).*besselh(n-g,k.*Rs).*am./4./mu./(1i).^(n-g-1)...
    .*dbesselj(n-g,k.*epsi)./dhankel(n-g,k.*epsi));
U(:,:,n) = Ui+Us;

end
Uf = sum(U,3);
%------------------------------
lp1 = subplot(1,1,1);
mq = 200;
Uf(isnan(z)) = NaN;
Uoc = Uf;
Uic = Uf;
tic = t;
xic = x;
xoc = x;
yic = y;
yoc = y;
Uoc(r < Ro) = NaN;
xoc(r < Ro) = NaN;
yoc(r < Ro) = NaN;
Uic(r >= Ro) = NaN;
ric = Ro.^2.*(Ri-epsi)./del.^2+...
    (Ro.^2+epsi.^2-2.*Ri.*Ro).*r./del.^2+(Ri-epsi).*r.^2./del.^2;
ric(r >= Ro) = NaN;
tic(r >= Ro) = NaN;
xic = ric.*cos(tic);
yic = ric.*sin(tic);
contourf(xoc,yoc,real(Uoc),mq,'edgecolor','none');
axis equal
hold on
colorbar
colormap hot 
hold on
contourf(xic,yic,real(Uic),mq,'edgecolor','none');
hold on
phit = 0:2*pi/1000:2*pi;    
plot(Ri.*cos(phit),Ri.*sin(phit),'k','LineWidth',1.0);
hold on
plot(Ro.*cos(phit),Ro.*sin(phit),'k','LineWidth',1.0);
%------------------------------
% [ir,it] = meshgrid(epsi:(Rs-epsi)/mi:Rs,-pi:2*pi/mi:pi);
% zi = ones(length(ir),length(it));
% ix = ir.*cos(it);
% iy = ir.*sin(it);
% [or,ot] = meshgrid(Rs:(Rt-Rs)/mo:Rt,-pi:2*pi/mo:pi);
% zo = ones(length(or),length(ot));
% ox = or.*cos(ot);
% oy = or.*sin(ot);
% zo((ox-Rs).^2+oy.^2<=delta.^2) = NaN;
% ox(isnan(zo)) = NaN;
% oy(isnan(zo)) = NaN;
%title(['$\frac{b}{B}=\,$',num2str(1/k),'$ \,.$'],'Interpreter','latex');
%ylabel('$\displaystyle \frac{\hat{\sigma}^{r\phi}_{\left(1\right)}}{\mu}$','rot',360,'Interpreter','latex')
%gr1=legend(lp1,'Inclusion','Location','northeast');
%set(gr1,'Interpreter','Latex');
% contourf(ox,oy,real(Uof),mq2,'edgecolor','none');
% Usi = (-1i).^(n-g).*besselh(n-g,k.*r).*exp(1i.*(n-g).*t)...
%     .*((-1).^(n-g).*besselh(n-g,k.*Rs).*am./4./mu./(1i).^(n-g-1)...
%     .*dbesselj(n-g,k.*epsi)./dhankel(n-g,k.*epsi));
% Uo(:,:,n) = Uio+Uso;
%------------------------------
