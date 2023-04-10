
# The Unified Cluster Catalog: towards a unique and comprehensive database of stellar clusters

<!-- MarkdownTOC -->

- Databases cross-match
  - New databases
  - Comments on new databases
    - New databases with no data
  - Old databases
    - HAO21
    - DIAS21
    - CG20
    - BICA19
    - MWSC
    - TARRICQ22
    - DIAS21+CG20+MWSC+TARRICQ22+HAO21+BICA19 matching
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


Clusters with some particularity:
* Melotte_111 has a very large radius


### Databases cross-match

We separate the analysis in two sections: *New databases* and *Old databases*.

The first section contains databases recently published (from 2019 up to
2023) which list newly discovered clusters that are not present in older
databases.

The second section is made up of databases each compiling a set of known
open clusters, all of them overlapping to some extent. The most recent one,
HAO21, is also the largest with more than 3000 entries. After cleaning
(particularly HAO21) and cross-matching these databases we are left with
XXXX unique open clusters.


#### New databases

Some of these DBs contain members data (i.e.: the stars selected as cluster
members) and all of them have parallax and proper motions center (useful for
our membership analysis)

* HUNT23   **TODO**
Improving the open cluster census. II. An all-sky cluster catalogue with Gaia
DR3, Hunt & Reffert (2021)
739 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2023arXiv230313424H/abstract), [Vizier]() ???

* QIN23
Hunting for Neighboring Open Clusters with Gaia DR3: 101 New Open Clusters
within 500 pc; Qin et al. (2023)
101 clusters (OSCN_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2023ApJS..265...12Q/abstract), Tables provided by Qin & Chen

* LI23  **REMOVED**
LISC Catalog of Star Clusters. II. High Galactic Latitude Open Clusters in
Gaia EDR3; Li & Mao (2023)
56 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2023ApJS..265....3L/abstract), [Zenodo](https://zenodo.org/record/7603419)

* LI22
LISC Catalog of Star Clusters. I. Galactic Disk Clusters in Gaia EDR3;
Li et al. (2022)
61 clusters (LISC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...19L/abstract), [Zenodo](https://zenodo.org/record/5705371#.YZPASbFdsrs)

* HE22_2
Unveiling hidden stellar aggregates in the Milky Way: 1656 new star clusters
found in Gaia EDR3; He et al. (2022)
1656 clusters (CWNU_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2022arXiv220908504H/abstract), [Vizier](https://cdsarc.cds.unistra.fr/ftp/vizier.submit/he22c/)

* HE22_1
A Blind All-sky Search for Star Clusters in Gaia EDR3: 886 Clusters within
1.2 kpc of the Sun; He et al. (2022)
270 clusters (CWNU_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2022ApJS..262....7H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/262/7)

* HE22
New Open-cluster Candidates Found in the Galactic Disk Using Gaia DR2/EDR3
Data, He et al. 2022
541 clusters (CWNU_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2022ApJS..260....8H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/260/8)

* HAO22
Newly detected open clusters in the Galactic disk using Gaia EDR3,
Hao et al. 2022
704 clusters (OC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2022A%26A...660A...4H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/660/A4)

* CASTRO22
Hunting for open clusters in Gaia EDR3: 628 new open clusters found with
OCfinder, Castro-Ginard et al. 2022
(*numbered from UBC 1001 in order to differentiate from the UBC clusters found
in Gaia DR2*)
628 clusters (UBC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.118C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/661/A118)

* HUNT21   **TODO**
Improving the open cluster census. I. Comparison of clustering algorithms
applied to Gaia DR2 data, Hunt & Reffert (2021)
41 clusters (PHOC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.104H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/646/A104)

* HE21
A catalogue of 74 new open clusters found in Gaia Data-Release 2,
He et al. (2021)
74 clusters (CWNU_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2021RAA....21...93H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/other/RAA/21.93)

* FERREIRA21 **TODO**
New star clusters discovered towards the Galactic bulge direction using Gaia
DR2, Ferreira et al. (2021)
34 clusters (UFMG_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502L..90F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/502/L90)

* CASTRO20
Hunting for open clusters in Gaia DR2: 582 new open clusters .. Galactic disc,
Castro-Ginard et al. 2020
570 clusters (UBC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A..45C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/635/A45)

* FERREIRA20
Discovery and astrophysical properties of Galactic open clusters in dense
stellar fields using Gaia DR2, Ferreira et al. (2020)
25 clusters (UFMG_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2021F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/496/2021)

* CASTRO19
Hunting for open clusters in Gaia DR2: the Galactic anticentre,
Castro-Ginard et al. 2019
53 clusters (UBC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2019A%26A...627A..35C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/627/A35)

* FERREIRA19 **TODO**
Three new Galactic star clusters discovered in the field of the open cluster
NGC 5999 with Gaia DR2, Ferreira et al. (2019)
3 clusters (UFMG_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.5508F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/483/5508)

* LIUPANG19
A Catalog of Newly Identified Star Clusters in Gaia DR2, Liu & Pang 2019
76 clusters (no acronym given, FOF/LP)
[Article](https://ui.adsabs.harvard.edu/abs/2019ApJS..245...32L/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/245/32)

* SIM19
207 New Open Star Clusters within 1 kpc from Gaia Data Release, Sim et al. 2019
207 clusters (UPK_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2019JKAS...52..145S/abstract), From Table 2 in online article
**138 of these clusters are included in Cantat-Gaudin & Anders (2020; CGA20) as
UPK since for 66 they were not able to compute memberships, and 3 were
duplicated UBC clusters. I thus remove these 3 clusters: UPK 19 (UBC 32),
UPK 172 (UBC 10a), and UPK 327 (UBC 88); note that CGA20 mistakenly wrote
UPK 176 instead of UPK 172**

* CASTRO18
A new method for unveiling open clusters in Gaia. New nearby open clusters
confirmed by DR2, Castro-Ginard et al. 2018
23 clusters (UBC_XX)
[Article](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..59C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/618/A59)


#### Comments on new databases

HE22_2 lists the cluster 1951 right on the field of the NGC 5139 (Omega
Centauri) globular cluster. No real cluster is apparent there


##### New databases with no data

* CHI23   **TODO** 
 **WARNING**: https://twitter.com/CantatGaudin/status/1638133660875456515
Blind Search of The Solar Neighborhood Galactic Disk within 5kpc: 1179 new
Star clusters found in Gaia DR3, Chi et al. (2023)
1179 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2023arXiv230310380C/abstract), [Vizier]() ???

* CHI23_2
Identifying 46 New Open Cluster Candidates in Gaia EDR3 Using a Hybrid
pyUPMASK and Random Forest Method; Chi et al. (2023)
46 clusters
[Article][https://ui.adsabs.harvard.edu/abs/2023ApJS..265...20C/abstract], email sent to weishoulin@astrolab.cn

* CHI23
LISC Catalog of Open Clusters.III. 83 Newly found Galactic disk open clusters
using Gaia EDR3; Chi et al. (2023)
83 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2023arXiv230208926C/abstract), email sent to zhongmuli@126.com



#### Old databases

* TARRICQ22
Structural parameters of 389 local open clusters, Tarricq et al. (2022)
467 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A..59T/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/659/A59)

* HAO21
Evolution of the local spiral structure of the Milky Way revealed by open
clusters, Hao et al. (2021)
3794 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2021A%26A...652A.102H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/652/A102)

* DIAS21
Updated parameters of 1743 open clusters based on Gaia DR2, Dias et al. 2021
1743 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2021MNRAS.504..356D), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/504/356)

* CG20
Painting a portrait of the Galactic disc with its stellar clusters,
Cantat-Gaudin et al. 2020
2017 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2020A%26A...640A...1C), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/640/A1/table1)

* BICA19
A Multi-band Catalog of 10978 Star Clusters ... in the Milky Way,
Bica et al. 2019
2913 clusters (OCs)
[Article](https://ui.adsabs.harvard.edu/abs/2019AJ....157...12B/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/AJ/157/12/table3)

* MWSC
Data retrieved from [HEASARC](https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html) selecting all clusters with 'class' equal to
"OPEN STAR CLUSTER"
2913 clusters (OPEN STAR CLUSTER)


##### HAO21

The original number of entries in HAO21 is 3794. I remove UBC, Liupang19,
He2020 (HE21 clusters), UFMG (FERREIRA20 clusters), 37 BH_XX clusters
present in CG20 (including BH_200 which is actually ESO_92_5),
36 vdBergh-Hagen_XXX clusters (either present in CG20 as BH_XX clusters or
duplicated in HAO21 as BH_XX clusters) and renamed the remaining 4 to BH_XX. I
also remove the remaining duplicated clusters and clusters listed three times.
The final number of clusters in HAO21 is 2771

Triple entries:
ESO_368_14 (x2) + ESO_368-14; Basel_11a (x2) + Basel_11A; ESO_368_11 (x2) +
ESO_368-11; ESO_130_06 (x2) + ESO_130-06

Duplicated clusters are assigned very different parameters in this work, which
warns us to be cautious about the values here presented.

HAO21 presents 1179 matches with CG20. Out of these:
- 301,171 (~25%, 15%) have plx  rel. differences >25%,>50%
- 65,48   (~6%, 4%)   have pmRA rel. differences >25%,>50%
- 87,58   (~7%, 5%)   have pmDE rel. differences >25%,>50%

Some clusters have wildly incorrect parameters. For example the cluster
Melotte_25 is assigned a parallax of 0.264 mas by HAO21 when its true value is
larger than 21 mas.

##### DIAS21

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

Changed all FSR clusters from "FSR 00XX" to "FSR_XX" and all ESO clusters from
"ESO 00XX-0X" to "ESO_XXXX_XX" to match the other DBs
Removed clusters LIUPANG19 (177; LP_XXX; DIAS21 lists clusters that are marked
in LIUPAG19 as "Class 3" or "need confirmation" but does not list some
confirmed clusters in LIUPAG19 such as LP_58), CASTRO-19-20 (287; UBC_XXX),
SIM19 (127, UPK_XXX), UMFG clusters present in FERREIRA20

##### CG20


Cantat-Gaudin et al. (2020) contains
2017 clusters (234129 stars with P>0.7); parameters for 1867.
Mean   = 125 stars per cluster (234129/1867)

1. What is the G limit used?
  The hard limit is 18, but there are a 15 stars up to ~19.6. These belong to
  Hyades/Melotte 25 (9) and Melotte 111 (6)
2. How were very large clusters processed?
  They were not processed with UPMASK. This is explained in Sect 2.1.

Changed all FSR clusters from "FSR 00XX" to "FSR_XX" and all ESO clusters from
"ESO 00XX-0X" to "ESO_XXXX_XX" to match the other DBs
Removed LIUPANG19 (35; LP_XXX), CASTRO-18-19-20 (550; UBC_XXX_),
SIM19 (138, UPK_XXX) clusters

##### BICA19

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

There are 18 UBC clusters in BICA19 with 9 of them apparently associated to
a previously identified cluster: Alessi 161|UBC 3, TRSG 4|UBC 1,
Teutsch 179|UBC 9, AT 4|Alessi-Teutsch 4|Alessi 26|TRSG 3|UBC 12,
Alessi 190|UBC 5, Alessi 145|UBC 13, Renou 23|UBC 8, Alessi 94|UBC 2,
Teutsch 182|UBC 6. We remove the remaining 9.

Fixed
"Carraro 1.MWSC 1829" --> "Carraro 1|MWSC 1829"
"Cernik 39" --> "Czernik 39"
de Wit 1 --> Wit 1 (to match MWSC)
JS 1 --> Juchert-Saloran 1

##### MWSC

Changed all FSR clusters from "FSR 00XX" to "FSR_XX" and all ESO clusters from
"ESO 00XX-0X" to "ESO_XXXX_XX" to match the other DBs

vdBergh-Hagen --> BH

##### TARRICQ22

Changed all FSR clusters from "FSR 00XX" to "FSR_XX" and all ESO clusters from
"ESO 00XX-0X" to "ESO_XXXX_XX" to match the other DBs.
Removed LIUPANG19 (35; LP_XXX), CASTRO-18-19-20 (550; UBC_XXX_),
SIM19 (138, UPK_XXX) clusters



##### DIAS21+CG20+MWSC+TARRICQ22+HAO21+BICA19 matching

After cross-matching these cleaned databases we have:

```
DIAS21 1134
CG20 1294
MWSC 2858
TARRICQ22 274
HAO21 2771
BICA19_all_names 2895

Clusters in None matched to clusters in DIAS21: 0
Clusters in None not matched to clusters in DIAS21: 0
Clusters in DIAS21 not matched to clusters in None: 1134

Clusters in DIAS21 matched to clusters in CG20: 1093
Clusters in DIAS21 not matched to clusters in CG20: 41
Clusters in CG20 not matched to clusters in DIAS21: 201
Clusters in cross-matched catalog so far: 1335

Clusters in DIAS21+CG20 matched to clusters in MWSC: 1119
Clusters in DIAS21+CG20 not matched to clusters in MWSC: 216
Clusters in MWSC not matched to clusters in DIAS21+CG20: 1739
Clusters in cross-matched catalog so far: 3074

Clusters in DIAS21+CG20+MWSC matched to clusters in TARRICQ22: 274
Clusters in DIAS21+CG20+MWSC not matched to clusters in TARRICQ22: 2800
Clusters in TARRICQ22 not matched to clusters in DIAS21+CG20+MWSC: 0
Clusters in cross-matched catalog so far: 3074

Clusters in DIAS21+CG20+MWSC+TARRICQ22 matched to clusters in HAO21: 2638
Clusters in DIAS21+CG20+MWSC+TARRICQ22 not matched to clusters in HAO21: 436
Clusters in HAO21 not matched to clusters in DIAS21+CG20+MWSC+TARRICQ22: 133
Clusters in cross-matched catalog so far: 3207

Clusters in DIAS21+CG20+MWSC+TARRICQ22+HAO21 matched to clusters in BICA19_all_names: 1931
Clusters in DIAS21+CG20+MWSC+TARRICQ22+HAO21 not matched to clusters in BICA19_all_names: 1276
Clusters in BICA19_all_names not matched to clusters in DIAS21+CG20+MWSC+TARRICQ22+HAO21: 998
Clusters in final cross-matched catalog: 4205
```

There are 999 clusters with neither parallax nor proper motions data and 1317
clusters with no parallax data. Most of these are the Ryu clusters.

There are 108 clusters with data only from HAO21, which is of low quality.

**Annotated clusters:**

- Loden 466 & Stock 14 appear to be the same cluster



### Cluster definition

Consider a point c in an n-dimensional space and an integer value m>0. We define
a 'cluster' as the collection of m elements with the smallest n-dimensional
euclidean distance to c.




### Annotated clusters

IC_2157:  two separated clusters?
IC_5146:  two separated clusters?
UBC 274 (https://ui.adsabs.harvard.edu/abs/2022A%26A...664A..31C/abstract)
NGC_1333: second binary cluster at ~(157.9, -21.3). Possible trail of stars
being ripped by IC_348
Mamajek_1: there appears to be another unidentified cluster at ~(300.4, -16)
NGC_2244: 2nd cluster in frame at ~(206.17, -2.3)?


### Clusters with issues

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




### Probably useful references

- [Gaia colour-magnitude diagrams of young open clusters.., Negueruela &
   de Burgos (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230300467N/abstract)
- [Discovery and description of two young open clusters in the primordial
   group of NGC 6871, Casado & Hendy (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp..104C/abstract)
- [Baumgardt catalog of globular clusters](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)

