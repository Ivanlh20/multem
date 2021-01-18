import numpy
import time
import multem
import multem.crystalline_materials

from math import log


def run():

    Z = 49
    charge = 0

    Rmin = 1e-02
    Rmax = 5.0
    nR = 512
    dlnR = log(Rmax / Rmin) / (nR - 1)
    R = Rmin * numpy.exp(numpy.arange(0, (nR - 1), 1) * dlnR)

    st = time.perf_counter()
    [f1, df1] = multem.vp(1, Z, charge, R)
    [f2, df2] = multem.vp(2, Z, charge, R)
    [f3, df3] = multem.vp(3, Z, charge, R)
    [f4, df4] = multem.vp(4, Z, charge, R)
    [f5, df5] = multem.vp(5, Z, charge, R)
    [f6, df6] = multem.vp(6, Z, charge, R)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # sub#plot(1, 2, 1)
    # hold on
    # #plot(R, f1, '-k', R, f2, '-b', R, f3, '-c', R, f4, '-m', R, f5, '-r', R, f6, '-g')
    # #set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title('Projected Atomic potential')
    # ylabel('$\displaystyle V(R)$','interpreter','latex','FontSize',14)
    # xlabel('$\mathbf{R}$','interpreter','latex','FontSize',12)
    # xlim([0 Rmax])
    # legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]')

    # sub#plot(1, 2, 2)
    # hold on
    # #plot(R, df1, '-k', R, df2, '-b', R, df3, '-c', R, df4, '-m', R, df5, '-r', R, df6, '-g')
    # #set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title('Derivative of the Atomic potential')
    # ylabel('$\displaystyle \frac{d V(R)}{dr}$','interpreter','latex','FontSize',14)
    # xlabel('$\mathbf{R}$','interpreter','latex','FontSize',12)
    # xlim([0 Rmax])
    # legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]')

    # #set(gcf,'units','normalized','outerposition',[0 0 1 1])


if __name__ == "__main__":
    run()
