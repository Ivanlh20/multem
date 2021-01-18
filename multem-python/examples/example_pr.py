import numpy
import multem
import multem.crystalline_materials

from math import log


def run():

    Z = 79
    occ = 1
    region = 0
    charge = 0

    rmin = 1e-05
    rmax = 0.1
    nr = 512
    dlnr = log(rmax / rmin) / (nr - 1)
    r = rmin * numpy.exp(numpy.arange(0, (nr - 1), 1) * dlnr)

    [f1, df1] = multem.pr(1, Z, charge, r)
    [f2, df2] = multem.pr(2, Z, charge, r)
    [f3, df3] = multem.pr(3, Z, charge, r)
    [f4, df4] = multem.pr(4, Z, charge, r)
    [f5, df5] = multem.pr(5, Z, charge, r)
    [f6, df6] = multem.pr(6, Z, charge, r)

    # sub#plot(1, 2, 1)
    # #plot(r, f1, '-k', r, f2, '-b', r, f3, '-c', r, f4, '-m', r, f5, '-r', r, f6, '-g')
    # #set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title('Electron density')
    # ylabel('$\displaystyle \rho(r)$','interpreter','latex','FontSize',14)
    # xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12)
    # xlim([0 rmax]) ylim([0 max(f6)])
    # legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]')

    # sub#plot(1, 2, 2)
    # #plot(r, df1, '-k', r, df2, '-b', r, df3, '-c', r, df4, '-m', r, df5, '-r', r, df6, '-g')
    # #set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title('Derivative of the Electron density')
    # ylabel('$\displaystyle \frac{d \rho(r)}{d r}$','interpreter','latex','FontSize',14)
    # xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12)
    # xlim([0 rmax]) ylim([min(df6) max(df6)])
    # legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]')

    # #set(gcf,'units','normalized','outerposition',[0 0 1 1])


if __name__ == "__main__":
    run()
