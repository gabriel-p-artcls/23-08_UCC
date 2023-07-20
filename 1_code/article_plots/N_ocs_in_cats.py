
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd


def main():
    """
    """
    years = (1771, 1900, 1995, 2015, 2022, 2023)
    values = (30, 650, 1200, 3000, 7000, 14000)
    # plt.bar(years, values, width=15)
    plt.plot(years, values, alpha=.75, marker='o', ms=10, color='maroon')
    plt.axvline(2016, ls=':', lw=2, c='k')
    # plt.text(1920, 5000, "Gaia release")

    plt.annotate(
        'Gaia release', xy=(2015, 5000), xytext=(1910, 5000), # fontsize=8,
        verticalalignment="center",
        # Custom arrow
        arrowprops=dict(arrowstyle='->', lw=.7))

    # plt.yscale('log')
    plt.xlabel("Year")#, fontsize=15)
    plt.ylabel("Catalogued OCs")#, fontsize=15)
    plt.savefig("catalogued_ocs.png", dpi=300)


if __name__ == '__main__':
    import scienceplots
    plt.style.use('science')
    main()
