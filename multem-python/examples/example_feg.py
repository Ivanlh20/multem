import numpy

import multem
import multem.crystalline_materials


def run():

    Z = 50
    occ = 1
    region = 0
    charge = 0

    gmin = 0
    gmax = 12
    ng = 512
    dg = (gmax - gmin) / (ng - 1)
    g = numpy.arange(gmin, gmax, dg)

    [f1, df1] = multem.feg(1, Z, charge, g)
    [f2, df2] = multem.feg(2, Z, charge, g)
    [f3, df3] = multem.feg(3, Z, charge, g)
    [f4, df4] = multem.feg(4, Z, charge, g)
    [f5, df5] = multem.feg(5, Z, charge, g)
    [f6, df6] = multem.feg(6, Z, charge, g)

    # sub#plot(1, 2, 1)
    ##plot(g, f1, "-k", g, f2, "-b", g, f3, "-c", g, f4, "-m", g, f5, "-r", g, f6, "-g")
    ## set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title("Electron Scattering factor")
    # ylabel("$\displaystyle f_e(g)$", "interpreter", "latex", "FontSize", 14)
    # xlabel("$\mathbf{g}$", "interpreter", "latex", "FontSize", 12)
    # xlim([0, gmax])
    # legend(
    #    "Doyle [0-4]",
    #    "Peng [0-4]",
    #    "Peng [0-12]",
    #    "Kirkland [0-12]",
    #    "Weickenmeier [0-12]",
    #    "Lobato [0-12]",
    # )

    # sub#plot(1, 2, 2)
    ##plot(
    #    g,
    #    df1,
    #    "-k",
    #    g,
    #    df2,
    #    "-b",
    #    g,
    #    df3,
    #    "-c",
    #    g,
    #    df4,
    #    "-m",
    #    g,
    #    df5,
    #    "-r",
    #    g,
    #    df6,
    #    "-g",
    # )
    ## set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title("Derivative of the electron Scattering factor")
    # ylabel(
    #    "$\displaystyle \frac{d f_e(g)}{dg}$", "interpreter", "latex", "FontSize", 14
    # )
    # xlabel("$\mathbf{g}$", "interpreter", "latex", "FontSize", 12)
    # xlim([0, gmax])
    # legend(
    #    "Doyle [0-4]",
    #    "Peng [0-4]",
    #    "Peng [0-12]",
    #    "Kirkland [0-12]",
    #    "Weickenmeier [0-12]",
    #    "Lobato [0-12]",
    # )

    # set(gcf,'units','normalized','outerposition',[0 0 1 1])


if __name__ == "__main__":
    run()
