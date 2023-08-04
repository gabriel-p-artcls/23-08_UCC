
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')


"""
"""
years = (1771, 1900, 1995, 2015, 2022, 2023)
values = (30, 650, 1200, 3000, 7000, 14000)
# plt.bar(years, values, width=15)
plt.plot(years, values, alpha=.75, marker='o', ms=10, color='maroon')
plt.axvline(2016, ls=':', lw=2, c='k')

plt.annotate(
    'Gaia release', xy=(2015, 5000), xytext=(1910, 5000), # fontsize=8,
    verticalalignment="center",
    # Custom arrow
    arrowprops=dict(arrowstyle='->', lw=.7))

# plt.yscale('log')
plt.xlim(1750, 2039)
plt.xlabel("Year")#, fontsize=15)
plt.ylabel("Catalogued OCs")#, fontsize=15)
plt.savefig("catalogued_ocs.png", dpi=300)
