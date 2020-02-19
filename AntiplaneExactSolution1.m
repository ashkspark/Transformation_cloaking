clear all; clc;
mu = 1.0;
lam = 1.0;
rho = 1.0;
omega = 30.28;
am = 20.0;
Rt = 1.0;
Rs = 0.7*Rt;
epsi = Rs/40;
mi = 100;
mo = mi;
delta = 0.07*Rt;
k = (rho*omega^2/mu)^(1/2);
s = 191;
g=(s+1)/2;
[ir,it] = meshgrid(epsi:(Rs-epsi)/mi:Rs,-pi:2*pi/mi:pi);
zi = ones(length(ir),length(it));
ix = ir.*cos(it);
iy = ir.*sin(it);
zi((ix-Rs).^2+iy.^2<=delta.^2) = NaN;
ix(isnan(zi)) = NaN;
iy(isnan(zi)) = NaN;
[or,ot] = meshgrid(Rs:(Rt-Rs)/mo:Rt,-pi:2*pi/mo:pi);
zo = ones(length(or),length(ot));
ox = or.*cos(ot);
oy = or.*sin(ot);
zo((ox-Rs).^2+oy.^2<=delta.^2) = NaN;
ox(isnan(zo)) = NaN;
oy(isnan(zo)) = NaN;
for n = 1:1:s
Uii = am./4./1i./mu.*besselh(n-g,k.*Rs).*besselj(n-g,k.*ir).*exp(1i.*(n-g).*it);
Uio = am./4./1i./mu.*besselh(n-g,k.*or).*besselj(n-g,k.*or).*exp(1i.*(n-g).*ot);
Usi = (-1i).^(n-g).*besselh(n-g,k.*ir).*exp(1i.*(n-g).*it)...
    .*((-1).^(n-g).*besselh(n-g,k.*Rs).*am./4./mu./(1i).^(n-g-1)...
    .*dbesselj(n-g,k.*epsi)./dhankel(n-g,k.*epsi));
Uso = (-1i).^(n-g).*besselh(n-g,k.*or).*exp(1i.*(n-g).*ot)...
    .*((-1).^(n-g).*besselh(n-g,k.*Rs).*am./4./mu./(1i).^(n-g-1)...
    .*dbesselj(n-g,k.*epsi)./dhankel(n-g,k.*epsi));
Ui(:,:,n) = Uii+Usi;
Uo(:,:,n) = Uio+Uso;
end
Uif = sum(Ui,3);
Uof = sum(Uo,3);
%------------------------------
lp1 = subplot(1,1,1);
mq1 = 100;
mq2 = mq1;
Uif(isnan(zi)) = NaN;
Uof(isnan(zo)) = NaN;
contourf(ix,iy,real(Uif),mq1,'edgecolor','none');
axis equal
hold on
contourf(ox,oy,real(Uof),mq2,'edgecolor','none');
colorbar
colormap hot 
hold on
phit = 0:2*pi/1000:2*pi;    
plot(epsi.*cos(phit),epsi.*sin(phit),'k','LineWidth',1.0);
%title(['$\frac{b}{B}=\,$',num2str(1/k),'$ \,.$'],'Interpreter','latex');
%ylabel('$\displaystyle \frac{\hat{\sigma}^{r\phi}_{\left(1\right)}}{\mu}$','rot',360,'Interpreter','latex')
%gr1=legend(lp1,'Inclusion','Location','northeast');
%set(gr1,'Interpreter','Latex');
%------------------------------
