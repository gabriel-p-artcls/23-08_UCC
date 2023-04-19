
# The Unified Cluster Catalog: towards a unique and comprehensive database of stellar clusters

<!-- MarkdownTOC -->

- Databases included in cross-match
  - KHARCHENKO12
  - CASTRO18
  - BICA19
  - SIM19
    - CASTRO19
  - LIUPANG19
  - FERREIRA19
  - CASTRO20
  - FERREIRA20
  - CANTAT20
  - HAO20
  - FERREIRA21
  - HE21
  - DIAS21
  - HUNT21
  - CASADO21
  - JAEHNIG21
  - SANTOS21
  - TARRICQ22
  - CASTRO22
  - HE22
  - HE22_1
  - HE22_2
  - HAO22
  - LI22
  - HUNT23
  - QIN23
  - LI23
  - CHI23_2
  - CHI23
  - CHI23_3
  - Small databases included in HUNT23
- Databases not included
  - HAO21
  - KOUNKEL20
- Cross-matching results
- Cluster definition
- Annotated clusters
- Clusters with issues
- Probably useful references

<!-- /MarkdownTOC -->

This project began on June 2022.


TODO:

- [x] Finish writing the fastMP code
  - [x] Add magnitude cut test <-- Not good results
  - [x] Add ABC classification
  - [x] Test t-snes dim reduction <-- does not help
- [ ] Add new databases [NASA/ADS Search](https://ui.adsabs.harvard.edu/search/filter_database_fq_database=AND&filter_database_fq_database=database%3A%22astronomy%22&filter_property_fq_property=AND&filter_property_fq_property=property%3A%22refereed%22&fq=%7B!type%3Daqp%20v%3D%24fq_database%7D&fq=%7B!type%3Daqp%20v%3D%24fq_property%7D&fq_database=(database%3A%22astronomy%22)&fq_property=(property%3A%22refereed%22)&q=abs%3A%22new%20open%20clusters%22&sort=date%20desc%2C%20bibcode%20desc&p_=0)
  - [ ] Download data for all databases including members
- [x] Generate a single database with all the necessary data
- [x] Make ucc.ar the landing site
- [x] Write script to fetch the cluster region + process with fastMP
- [x] Write script to generate plots
- [x] Write a Jupyter notebook template
- [ ] Write the landing page for the site
- [x] Write the template site for the clusters (https://app.aavso.org/vsp/ ?)
- [ ] Store the un-edited DBs to retrieve the fundamental parameters
  - [ ] Write script to retrieve the fundamental parameters
- [ ] Write script to name each cluster
- [ ] Write script to generate each template page + plots + notebook
- [ ] Write the article



## Databases included in cross-match

Before combining the databases we must first sanitize them so that they
can be read appropriately. These are the databases and the changes made to each:

### KHARCHENKO12
Global survey of star clusters in the Milky Way. I. The
Pipeline and ...; [Kharchenko et al. 2012](https://ui.adsabs.harvard.edu/abs/2012A%26A...543A.156K), [HEASARC](https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html)

I retrieved the data from the HEASARC service selecting all clusters with
'class' equal to "OPEN STAR CLUSTER", resulting in 2858 entries.
This DB lists the FSR clusters as "FSR XXXX" with leading zeros and the ESO
clusters as "ESO XXX-XX" with leading zeroes.

vdBergh-Hagen --> VDBH per CDS recommendation

### CASTRO18
A new method for unveiling open clusters in Gaia. New nearby open
clusters confirmed by DR2, [Castro-Ginard et al. 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..59C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/618/A59)

The table contains 23 UBC clusters.

### BICA19
A Multi-band Catalog of 10978 Star Clusters ... in the Milky Way;
2913 clusters (OCs); [Bica et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....157...12B/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/AJ/157/12/table3)

The Vizier table contain 10978 entries. We keep only those with Class1 OC (open
cluster) or OCC (open cluster candidate); this reduces the list to 3564 entries.
This DB lists the FSR clusters as "FSR XXX" and the ESO clusters as
"ESO XXX-XX" with no leading zeroes. I added 'RA_ICRS' and
'DE_ICRS' columns to math the rest of the DBs.

The list contains 19 duplicates and one entry where the same name is listed
twice (MWSC 2776). 
- For the single entry with the same name listed twice, we remove one of the
  repeated names.
- For (FSR 523, FSR 847, FSR 436) we remove the duplicates that showed as
  single entries
- ESO 393-3: removed name from both entries (cluster not found in CDS)
- MWSC 1025, 1482, 948, 3123, 1997, 1840, 442, 1808, 2204: removed name from
  both entries (clusters not found in KHARCHENKO12)
- ESO 97-2: removed from Loden 848 as it matches the position of Loden 894
  according to CDS
- FSR 972, OCL 344, Collinder 384, FSR 179: removed from both entries, it does
  not show in CDS or anywhere else
- MWSC 206: removed entry that also showed FSR 60 as the coordinates for FSR 60
  are a better match in KHARCHENKO12 for the entry with the single FSR 60 name

Fixed:
FSR 429.MWSC 3667 --> FSR 429,MWSC 3667
Carraro 1.MWSC 1829 --> Carraro 1,MWSC 1829
Cernik 39 --> "Czernik 39
FSR343 --> FSR 343
ESO456-13 --> ESO 456-13
de Wit 1 --> Wit 1 (to match KHARCHENKO12)
JS 1 --> Juchert-Saloran 1 (to match KHARCHENKO12)
ESO 589-26,MW --> ESO 589-26
Messineo 1,Cl 1813-18,SAI 126, --> Remove comma at the end
BH --> VDBH per CDS recommendation

This is the only DB that lists the Ryu & Lee (2018) clusters. The original
article claims to have found 721 new OCs (923 minus 202 embedded). BICA19 (page
11) says that the Ryu & Lee article lists 719 OCs (921 minus 202 embedded).
BICA19 lists in its Vizier table only 711 Ryu OCs, 4 of which are listed with
alternative names (Teutsch J1814.6-2814|Ryu 563, Quartet|Ryu 858,
GLIMPSE 70|Mercer 70|Ryu 273, LS 468|La Serena 468|Ryu 094). Hence there are
707 Ryu clusters in the final BICA19 Vizier table.

For the Ryu clusters:
*The average angular radius of the clusters is 1'.31±0 60.
More specifically, 902 (98%) clusters are smaller than 3′, and
823 (89%) clusters are even smaller than 2′.* [Ryu & Lee (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...856..152R/abstract)

### SIM19
207 New Open Star Clusters within 1 kpc from Gaia Data Release,
[Sim et al. 2019](https://ui.adsabs.harvard.edu/abs/2019JKAS...52..145S/abstract), Data taken from Table 2 in online article

The table contains 207 UPK clusters. Three of these clusters were identified
as duplicated UBC clusters IN CG20. I thus remove these 3 clusters: UPK 19
(UBC 32), UPK 172 (UBC 10a), and UPK 327 (UBC 88); note that CGA20 mistakenly
wrote UPK 176 instead of UPK 172. Added (ra, dec) columns and a plx column
estimated as 1000/dist_pc

#### CASTRO19
Hunting for open clusters in Gaia DR2: the Galactic anticentre,
[Castro-Ginard et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...627A..35C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/627/A35)

The table contains 53 UBC clusters.

### LIUPANG19
A Catalog of Newly Identified Star Clusters in Gaia DR2,
[Liu & Pang 2019](https://ui.adsabs.harvard.edu/abs/2019ApJS..245...32L/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/245/32)

The table contains 76 clusters with no acronym given. I added 'FoF_' to match
HUNT23.

### FERREIRA19
Three new Galactic star clusters discovered in the field of
the open cluster NGC 5999 ..., [Ferreira et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.5508F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/483/5508)

The table contains 3 UFMG clusters. Added (ra, de) columns.

### CASTRO20
Hunting for open clusters in Gaia DR2: 582 new open clusters ..
Galactic disc, [Castro-Ginard et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A..45C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/635/A45)

The table contains 570 UBC clusters.

### FERREIRA20
Discovery and astrophysical properties of Galactic open
clusters in dense..., [Ferreira et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2021F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/496/2021)

The table contains 25 UFMG clusters. Added (ra, de) columns.

### CANTAT20
Painting a portrait of the Galactic disc with its stellar clusters;
[Cantat-Gaudin et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...640A...1C), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/640/A1/table1)

The table contains 2017 clusters. This DB lists the FSR clusters as "FSR_XXXX"
with leading zeros and the ESO clusters as "ESO_XXX_XX" with leading zeroes.

BH --> VDBH per CDS recommendation
LP_ --> FoF_ to match the original work LIUPANG19

### HAO20
Sixteen Open Clusters Discovered with Sample-based Clustering Search of Gaia
DR2, [Hao et al. 2020](https://ui.adsabs.harvard.edu/abs/2020PASP..132c4502H/abstract), Data from Table 2

This table lists 16 clusters with no acronym. Used 'HXWHB_' to match HUNT23

### FERREIRA21
New star clusters discovered towards the Galactic bulge
direction using Gaia DR2, [Ferreira et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502L..90F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/502/L90)

The table contains 34 UFMG clusters.

### HE21
A catalogue of 74 new open clusters found in Gaia Data-Release 2,
[He et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021RAA....21...93H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/other/RAA/21.93)

The table contains 74 with no acronym, added 'CWNU_'. Added (ra, de) columns.

### DIAS21
Updated parameters of 1743 open clusters based on Gaia DR2,
[Dias et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.504..356D), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/504/356)

The table contains 1743 clusters. This DB lists the FSR clusters as "FSR_XXXX"
with leading zeros and the ESO clusters as "ESO_XXX_XX" with leading zeroes.

The table lists 177 LIUPANG19 clusters because it includes clusters not
listed as new by the authors. We remove all except those listed as new in
LIUPANG19. Cluster 866 was duplicated so one entry was removed. Changed
'LP_' to 'FoF_' for consistency.

BH --> VDBH per CDS recommendation

This DB contains wrong RA data for several clusters (differences with HAO21:
`10>1 deg`, `14>0.5 deg`, `86>0.1 deg`)

The `r50` column is listed with 'pc' units in Vizier, but it is degrees unit.
It also shows 4 clusters with very large `r50` values which could be listed
in pc units:

```
Cluster      DIAS21  CG20   BICA19
----------------------------------
Berkeley_58  32.969  0.06   0.058
Blanco_1     13.218  0.699  0.833
NGC_7789     9.324   0.211  0.133
Berkeley_59  3.097   0.137  nan
```

### HUNT21
Improving the open cluster census. I. Comparison of clustering algorithms
applied to Gaia DR2 data, [Hunt & Reffert (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.104H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/646/A104)

This table lists 41 'PHOC_' clusters.

### CASADO21
New open clusters found by manual mining of data based on Gaia DR2,
[Casado (2021)](https://ui.adsabs.harvard.edu/abs/2021RAA....21..117C/abstract), Data from Table 1 in article

This table lists 20 'Casado_' clusters.

### JAEHNIG21
Membership Lists for 431 Open Clusters in Gaia DR2 Using Extreme Deconvolution
Gaussian Mixture Models, [Jaehnig et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...923..129J/abstract), Data from Table 1 in article

This table lists 11 'XDOCC-' clusters. Changed 'XDOCC-0Y' to 'XDOCC-Y' to match
HUNT23.

### SANTOS21
Canis Major OB1 stellar group contents revealed by Gaia,
[Santos-Silva et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.1033S/abstract), Data from Table 1 in article

This table lists 5 'CMa-' clusters.

### TARRICQ22
Structural parameters of 389 local open clusters,
[Tarricq et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A..59T/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/659/A59)

The table contains 467 clusters. This DB lists the FSR clusters as "FSR_XXX"
with leading zeros and the ESO clusters as "ESO_XXX_XX" with leading zeroes.

BH_99 --> VDBH_99
LP_ --> FoF_

### CASTRO22
Hunting for open clusters in Gaia EDR3: 628 new open clusters
found with OCfinder, [Castro-Ginard et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.118C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/661/A118)

The table contains 628 UBC clusters numbered from UBC 1001 in order to
differentiate them from the UBC clusters found using Gaia DR2 in the previous
articles.

### HE22
New Open-cluster Candidates Found in the Galactic Disk Using Gaia
DR2/EDR3 Data, [He et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJS..260....8H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/260/8)

The table contains 541 with no acronym, added 'CWNU_'.

### HE22_1
A Blind All-sky Search for Star Clusters in Gaia EDR3: 886
Clusters within 1.2 kpc of the Sun, [He et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJS..262....7H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/262/7)

The table contains 270 'CWNU_' clusters.

### HE22_2
Unveiling hidden stellar aggregates in the Milky Way: 1656 new star
clusters found in Gaia EDR3, [He et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220908504H/abstract), [Vizier](https://cdsarc.cds.unistra.fr/ftp/vizier.submit/he22c/)

The table contains 1656 with no acronym, added 'CWNU_'. Added (ra, de) columns.
This DB lists the cluster 1951 right on the field of the NGC 5139 (Omega
Centauri) globular cluster. No real cluster is apparent there

### HAO22
Newly detected open clusters in the Galactic disk using Gaia EDR3,
[Hao et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...660A...4H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/660/A4)

This table lists 704 'OC_' clusters.

### LI22
LISC Catalog of Star Clusters. I. Galactic Disk Clusters in Gaia EDR3,
[Li et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...19L/abstract), [Zenodo](https://zenodo.org/record/5705371#.YZPASbFdsrs)

This table lists 61 'LISC_' clusters.

### HUNT23
Improving the open cluster census. II. An all-sky cluster catalogue
with Gaia DR3, [Hunt & Reffert (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230313424H/abstract), [Vizier]() ???

The table contains 7167 clusters. This DB lists the FSR clusters as "FSR_XXX"
with no leading zeros and the ESO clusters as "ESO_XXX-XX" with no leading
zeroes.

Added the names of the 'HSC_' clusters to the 'name_all' column that were
missing, except for line 3400 (HSC_1995) were a cluster was already identified.
Added the name of the NGC_2533 cluster in the "name_all" column (4587 line)
which was also missing.

Out of the 563 'Theia' moving groups listed we keep the 252 with no other
objects listed nearby.

Fixed:
ESO 589-26,MW --> ESO 589-26 (dragged from BICA19)
ESO_429-429 --> ESO_429-02 (according to CDS coords)
BH --> VDBH per CDS recommendation

### QIN23
Hunting for Neighboring Open Clusters with Gaia DR3: 101 New Open
Clusters within 500 pc, [Qin et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJS..265...12Q/abstract), Tables provided by Qin & Chen

This table lists 101 'OSCN_' clusters.

### LI23
LISC Catalog of Star Clusters. II. High Galactic Latitude Open Clusters in
Gaia EDR3; [Li & Mao (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJS..265....3L/abstract), [Zenodo](https://zenodo.org/record/7603419)

This table lists 56 'LISC' clusters but only 35 are kept as "real" objects.
The parallax distances are in very bad agreement with the the estimated
distance moduli. These appear to either be MC clusters or not real clusters
at all. HUNT23 recovers 0% of these clusters.

### CHI23_2
Identifying 46 New Open Cluster Candidates in Gaia EDR3 Using a Hybrid
pyUPMASK and Random Forest Method; [Chi et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJS..265...20C/abstract), data taken
from Table in the article (IOP)

This table lists 46 clusters with no acronym, 'CWWL_' was added to match
HUNT23.

### CHI23
LISC Catalog of Open Clusters.III. 83 Newly found Galactic disk open clusters
using Gaia EDR3; [Chi et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230208926C/abstract)

The article mentions 83 clusters but only 82 are visible in the PDF. I sent an
email to zhongmuli@126.com but never got an answer. Added 'LISC-III' to the
names to match HUNT23. Added (ra, dec) columns.

### CHI23_3
**WARNING**: https://twitter.com/CantatGaudin/status/1638133660875456515
Blind Search of The Solar Neighborhood Galactic Disk within 5kpc: 1179 new
Star clusters found in Gaia DR3, [Chi et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230310380C/abstrac), Table requested to
authors

The table lists 1179 clusters. Added CWWDL following the convention by HUNT23
for DBs with no acronyms.


### Small databases included in HUNT23

These 6 articles are small completely covered by HUNT23

Qin et al. (2021), Casado & Hendy (2023), Anders et al. (2022),
Bastian (2019), Tian (2020), Zari et al. (2018)


## Databases not included

### HAO21
Evolution of the local spiral structure of the Milky Way revealed by
open clusters, [Hao et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...652A.102H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/652/A102)

The original number of entries in HAO21 is 3794. This database contains many
errors. It contains dozens of duplicated entries and even some that are listed
thrice, e.g: ESO_368_14 (x2) + ESO_368-14; Basel_11a (x2) + Basel_11A;
ESO_368_11 (x2) + ESO_368-11; ESO_130_06 (x2) + ESO_130-06. Furtheremore
duplicated clusters are assigned very different fundamental parameters, which
warns us to be cautious about the values here presented.

HAO21 presents 1179 matches with CG20. Out of these:
- 301,171 (~25%, 15%) have plx  rel. differences >25%,>50%
- 65,48   (~6%, 4%)   have pmRA rel. differences >25%,>50%
- 87,58   (~7%, 5%)   have pmDE rel. differences >25%,>50%

Some clusters have wildly incorrect parameters. For example the cluster
Melotte_25 is assigned a parallax of 0.264 mas by HAO21 when its true value is
larger than 21 mas.

For these reasons we exclude this database from the list.

### KOUNKEL20
Untangling the Galaxy. II. Structure within 3 kpc,
[Kounkel et al. 2022](https://ui.adsabs.harvard.edu/abs/2020AJ....160..279K/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AJ/160/279)

This table lists moving groups rather than open clusters. Kounkel private
communication:
"*All of these populations are moving groups. The word "cluster" is so nebulous
nowadays, we don't have a proper definition of what it means any more, outside
of those classicaly identified.*"

We thus include in our catalogue only a fraction of those listed by HUNT23.



## Cross-matching results

After cross-matching the sanitized databases we have:

```
1 KHARCHENKO12 2858
2 CASTRO18 23
3 BICA19 3560
4 CASTRO19 53
5 SIM19 204
6 LIUPANG19 76
7 FERREIRA19 3
8 CANTAT20 2017
9 CASTRO20 570
10 FERREIRA20 25
11 HAO20 16
12 DIAS21 1604
13 CASADO21 20
14 FERREIRA21 34
15 HUNT21 41
16 JAEHNIG21 11
17 SANTOS21 5
18 HE21 74
19 CASTRO22 628
20 TARRICQ22 467
21 LI22 61
22 HE22 541
23 HE22_1 270
24 HE22_2 1656
25 HAO22 704
26 HUNT23 6856
27 QIN23 101
28 LI23 35
29 CHI23_2 46
30 CHI23 82
31 CHI23_3 1179

23820 clusters in all DBs
Clusters in final cross-matched catalog: 13967
```



## Cluster definition

Consider a point c in an n-dimensional space and an integer value m>0. We define
a 'cluster' as the collection of m elements with the smallest n-dimensional
euclidean distance to c.




## Annotated clusters

IC_2157:  two separated clusters?
IC_5146:  two separated clusters?
UBC 274 (https://ui.adsabs.harvard.edu/abs/2022A%26A...664A..31C/abstract)
NGC_1333: second binary cluster at ~(157.9, -21.3). Possible trail of stars
being ripped by IC_348
Mamajek_1: there appears to be another unidentified cluster at ~(300.4, -16)
NGC_2244: 2nd cluster in frame at ~(206.17, -2.3)?
Loden 466 & Stock 14 appear to be the same cluster
Melotte_111 has a very large radius


## Clusters with issues

1. BDSB93: CMD (not too) contaminated bc of too many members

2. Collinder_419: ~200 extra members assigned by fastMP

3. Mamajek_1: too many members: The cluster is close but very concentrated.
Reducing the size of the frame from 20 to 10 fixes the issue

4. NGC_2645: too many members. Pismis 8 in frame but removed. There appears to be
a portion of the clusters separated in coordinates, centered at ~(264.9, -2.88)
I estimate the number of members at around 120. A larger C_thresh (x5) provides
the correct number of members

5. Berkeley_82: (not) too many members. Manually I estimate ~110, the code
gives ~180

6. FSR_0158: too many members. Manually I estimate ~170, the code gives ~380.
A larger C_thresh (x2) provides the correct number of members

7. UPK_542: moves + too few members. Fixing the xy+vpd centers fixes the issue

8. UBC_610: too few members. Fixing the xy+vpd centers fixes the issue

9. UBC_591: moves + too few members. Manually --> ~50 members. fastMP also
estimates ~50 members (Ripley). Fixing the xy+vpd centers and selecting
a smaller prob_min (to match the number of members estimated) fixes the issue

10. UBC_493: moves to another cluster. Fixing  xy+vpd coords works but almost
400 members are estimated when manually I estimate ~100

11. UBC_55: moves to another cluster. Fixing  xy+vpd coords works
UBC_381: moves to another cluster. Fixing  xy+vpd coords works
Czernik_20: moves to another cluster. Fixing  xy+vpd coords works
UBC_549: moves to another cluster. Fixing  xy+vpd coords works
DBSB_100: moves to another cluster. Fixing  xy+vpd coords works
Gulliver_56: moves plx center. Fixing  xy+vpd coords works

12. ASCC_85: moves to another cluster. Fixing  xy+vpd coords works but estimates
~20 members when manually I find ~90

13. UBC_100: moves to another cluster. The clusters i very close to the GC and
contains ~1000000 stars Fixing the plx works

14. UBC_663: moves to another cluster. Fixing  xy+vpd+plx coords works
UBC_277: moves to another cluster

15. DBSB_60: moves to another cluster. Fixing xy+vpd fixes the position but gives
too many members. No cluster there?

16. Pismis_19: traces of NGC_5617. Clusters are too close to one another
Gulliver_29: contaminated by another cluster? Apparently not
NGC_3324: weird CMD. Apparently right
UBC_323: 2nd cluster in frame? Apparently not
UPK_621: CMD different from CG20. Fine apparently
UBC_603: moves plx center? Apparently not (wrong in CG)
NGC_2183: more clusters in frame. Apparently not. Cluster disrupted?




## Probably useful references

- [Gaia colour-magnitude diagrams of young open clusters.., Negueruela &
   de Burgos (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230300467N/abstract)
- [Discovery and description of two young open clusters in the primordial
   group of NGC 6871, Casado & Hendy (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp..104C/abstract)
- [Baumgardt catalog of globular clusters](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)

