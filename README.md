
# The Unified Cluster Catalog: towards a unique and comprehensive database of stellar clusters

<!-- MarkdownTOC -->

- converting from .py to .ipynb
- How can I run notebooks of a Github project in Google Colab?
  - Databases cross-match
    - New databases
      - New databases with no data
    - Old databases
      - HAO21
      - DIAS21
      - CG20
      - BICA19
      - MWSC
      - MWSC+DIAS21+CG20+BICA19 matching
  - Cluster definition
  - IAU naming convention
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
- [ ] Generate a single database with all the necessary data
- [x] Write script to fetch the cluster region + process with fastMP
- [ ] Write script to generate plots
- [x] Make ucc.ar the landing site
- [ ] Write the landing page for the site
- [ ] Write the template site for the clusters (https://app.aavso.org/vsp/)
- [ ] Write a Jupyter notebook for each cluster (open in Colab)?
- [ ] Write the article


Clusters with some particularity:
* Melotte_111 has a very large radius



# converting from .py to .ipynb

https://stackoverflow.com/questions/62510114/converting-from-py-to-ipynb
https://stackoverflow.com/questions/23292242/converting-to-not-from-ipython-notebook-format
jupytext: https://jupytext.readthedocs.io/en/latest/using-cli.html

# How can I run notebooks of a Github project in Google Colab?

https://stackoverflow.com/questions/62596466/how-can-i-run-notebooks-of-a-github-project-in-google-colab

* Change the domain from 'github.com' to 'githubtocolab.com
* Use:
https://colab.research.google.com/github/{user}/{repo-name}/blob/master/{path}/{notebook-name}.ipynb



### Databases cross-match

We separate the analysis in two sections: *New databases* and *Old databases*.

The first section contains databases recently published (from 2019 up to
2023) which list newly discovered clusters that are not present in older
databases.

The second section is made up of databases each compiling a set of known
open clusters, all of them overlapping to some extent. The most recent one,
HAO21, is also the largest with more than 3000 entries. After cleaning
(particularly HAO21) and cross-matching these 4 databases we are left with
XXXX unique open clusters.


#### New databases

Some of these DBs contain members data (i.e.: the stars selected as cluster
members) and all of them have parallax and proper motions center (useful for
our membership analysis)

* QIN23
Hunting for Neighboring Open Clusters with Gaia DR3: 101 New Open Clusters
within 500 pc; Qin et al. (2023)
101 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2023ApJS..265...12Q/abstract), Tables provided by Qin & Chen

* LI23  **REMOVE??**
LISC Catalog of Star Clusters. II. High Galactic Latitude Open Clusters in
Gaia EDR3; Li & Mao (2023)
56 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2023ApJS..265....3L/abstract), [Zenodo](https://zenodo.org/record/7603419)

* LI22
LISC Catalog of Star Clusters. I. Galactic Disk Clusters in Gaia EDR3;
Li et al. (2022)
61 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...19L/abstract), [Zenodo](https://zenodo.org/record/5705371#.YZPASbFdsrs)

* HE22_2
Unveiling hidden stellar aggregates in the Milky Way: 1656 new star clusters
found in Gaia EDR3; He et al. (2022)
1656 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022arXiv220908504H/abstract), [Vizier](https://cdsarc.cds.unistra.fr/ftp/vizier.submit/he22c/)

* HE22_1
A Blind All-sky Search for Star Clusters in Gaia EDR3: 886 Clusters within
1.2 kpc of the Sun; He et al. (2022)
270 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022ApJS..262....7H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/262/7)

* HE22
New Open-cluster Candidates Found in the Galactic Disk Using Gaia DR2/EDR3
Data, He et al. 2022
541 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022ApJS..260....8H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/260/8)

* HAO22
Newly detected open clusters in the Galactic disk using Gaia EDR3,
Hao et al. 2022
704 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022A%26A...660A...4H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/660/A4)

* CASTRO22
Hunting for open clusters in Gaia EDR3: 628 new open clusters found with
OCfinder, Castro-Ginard et al. 2022
628 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.118C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/661/A118)

* HE21
A catalogue of 74 new open clusters found in Gaia Data-Release 2,
He et al. (2021)
74 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2021RAA....21...93H/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/other/RAA/21.93)

* CASTRO20
Hunting for open clusters in Gaia DR2: 582 new open clusters .. Galactic disc,
Castro-Ginard et al. 2020
570 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A..45C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/635/A45)

* FERREIRA20
Discovery and astrophysical properties of Galactic open clusters in dense
stellar fields using Gaia DR2, Ferreira et al. (2020)
25 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2021F/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/496/2021)

* CASTRO19
Hunting for open clusters in Gaia DR2: the Galactic anticentre,
Castro-Ginard et al. 2019
53 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2019A%26A...627A..35C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/627/A35)

* LIUPANG19
A Catalog of Newly Identified Star Clusters in Gaia DR2, Liu & Pang 2019
76 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2019ApJS..245...32L/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/245/32)

* SIM19
207 New Open Star Clusters within 1 kpc from Gaia Data Release, Sim et al. 2019
207 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2019JKAS...52..145S/abstract), From Table 2 in online article
**138 of these clusters are included in Cantat-Gaudin & Anders (2020; CGA20) as
UPK since for 66 they were not able to compute memberships, and 3 were
duplicated UBC clusters. I thus remove these 3 clusters: UPK 19 (UBC 32),
UPK 172 (UBC 10a), and UPK 327 (UBC 88); note that CGA20 mistakenly wrote
UPK 176 instead of UPK 172**

* CASTRO18
A new method for unveiling open clusters in Gaia. New nearby open clusters
confirmed by DR2, Castro-Ginard et al. 2018
23 clusters
[Article](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..59C/abstract), [Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/618/A59)


##### New databases with no data

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

This DB contains duplicated entries for 116 clusters and 3 clusters that are
listed three times. We kept the clusters with the largest number of members,
checking that their proper motions were consistent with those from CG20, or
those from DIAS21 if CG20 proper motions were not available. In case proper
motions are very discrepant, we removed both clusters. If no proper motions
are available, we selected the one with the largest number of members.

Duplicated clusters like 'vdBergh-Hagen_34' have very discrepant values
throughout the literature or with other catalogues and were removed from
this catalog.

The original number of entries in HAO21 is 3794, and the final number of
entries is 3666 (`3794-3666 = 128 = 116+3*2+6`).

Removed triple entries:
```
2721,Basel_11a,109.281,0.033,-13.97,0.027,228.262,-0.773,0.713,0.03,1,-1.56,0.048,0.918,0.078,7.8,0.179,       ,      ,     ,  ,27,Dias et al. 2002         ,Basel 11a
686,Basel_11A,109.3,0.075,-13.995,0.062,228.292,-0.768,0.71,0.032,1,-1.597,0.016,0.974,0.019,7.42,      ,       ,      ,     ,  ,35,MWSC                     ,Basel 11A
1129,vdBergh-Hagen_87,151.173,0.097,-55.408,0.031,280.749,0.112,0.417,0.023,1,-6.176,0.071,2.692,0.05,7.3,      ,       ,      ,     ,  ,31,MWSC                     ,Cl VDBH 87
2658,BH_87,151.09,0.038,-55.423,0.018,280.719,0.072,0.456,0.026,1,-6.017,0.224,3.134,0.177,8.4,      ,       ,      ,     ,  ,30,Dias et al. 2002         ,Cl VDBH 87
1633,vdBergh-Hagen_217,259.08,0.044,-40.803,0.041,346.795,-1.508,0.396,0.025,1,-0.638,0.081,-1.819,0.077,7.7,      ,       ,      ,     ,  ,42,MWSC                     ,Cl VDBH 217
2701,BH_217,259.072,0.026,-40.827,0.026,346.772,-1.517,0.352,0.029,1,-0.705,0.017,-1.881,0.03,7.7,      ,       ,      ,     ,  ,45,Dias et al. 2002         ,Cl VDBH 217
```
Removed clusters:
```
943,vdBergh-Hagen_34,127.832,0.1,-44.482,0.074,262.569,-2.947,0.555,0.026,1,-4.148,0.275,4.494,0.209,8.48,      ,       ,      ,     ,  ,32,MWSC                     ,Cl VDBH 34
2817,BH_34,127.908,0.155,-44.528,0.101,262.638,-2.931,0.277,0.02,1,-3.211,0.054,3.982,0.057,8.88,0.204,45.45,0.81,0.1,3,64,Dias et al. 2002         ,Cl VDBH 34
1186,vdBergh-Hagen_99,159.454,0.332,-59.194,0.126,286.553,-0.639,0.331,0.015,1,-6.036,0.148,2.803,0.09,8,      ,       ,      ,     ,  ,87,MWSC                     ,Cl VDBH 99
2808,BH_99,159.443,0.174,-59.21,0.089,286.556,-0.656,0.355,0.02,1,-6.254,0.048,2.842,0.026,7.605,      ,-0.19,1.01,0.07,9,214,Dias et al. 2002         ,Cl VDBH 99
2867,BH_111,167.29,0.036,-63.836,0.015,291.944,-3.168,0.337,0.024,1,-6.863,0.253,2.301,0.115,8.35,      ,       ,      ,     ,  ,34,Dias et al. 2002         ,Cl VDBH 111
1227,vdBergh-Hagen_111,167.254,0.217,-63.842,0.103,291.931,-3.18,0.311,0.021,1,-6.13,0.102,2.031,0.068,8.35,      ,       ,      ,     ,  ,105,MWSC                     ,Cl VDBH 111
```
The fact that the same clusters are assigned wildly different parameters in this
works warns us to be cautious about the values here presented.

HAO21 presents 1807 matches with CG20. Out of these:
- 99,164  (~5.5%,9%)  have pmRA rel. differences >25%,>10%
- 118,179 (~6.5%,10%) have pmDE rel. differences >25%,>10%
- 383,662 (~21%,37%)  have plx rel. differences >25%,>10%

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

##### CG20

* Cantat-Gaudin & Anders (2020)

1481 clusters (435833 stars with P>0.01); G_max=18
Mean   = 294 stars per cluster (435833/1481)

* Cantat-Gaudin et al. (2020)

2017 clusters (234129 stars with P>0.7); parameters for 1867.
Mean   = 125 stars per cluster (234129/1867)

1. What is the G limit used?
  The hard limit is 18, but there are a 15 stars up to ~19.6. These belong to
  Hyades/Melotte 25 (9) and Melotte 111 (6)
2. How were very large clusters processed?
  They were not processed with UPMASK. This is explained in Sect 2.1.

Changed all FSR clusters from "FSR 00XX" to "FSR_XX" and all ESO clusters from
"ESO 00XX-0X" to "ESO_XXXX_XX" to match the other DBs

##### BICA19

This is the only DB that lists the Ryu & Lee (2018) clusters. The original
article claims to have found 721 new OCs (923 minus 202 embedded). BICA19 (page
11) says that the Ryu & Lee article lists 719 OCs (921 minus 202 embedded).
BICA19 lists in its Vizier table only 711 Ryu OCs, 4 of which are listed with
alternative names (Teutsch J1814.6-2814|Ryu 563, Quartet|Ryu 858,
GLIMPSE 70|Mercer 70|Ryu 273, LS 468|La Serena 468|Ryu 094). Hence there are
707 Ryu clusters in the final BICA19 Vizier table.

Fixed "Carraro 1.MWSC 1829" --> "Carraro 1|MWSC 1829"
Fixed "Cernik 39" --> "Czernik 39"

##### MWSC

Changed all FSR clusters from "FSR 00XX" to "FSR_XX" and all ESO clusters from
"ESO 00XX-0X" to "ESO_XXXX_XX" to match the other DBs

vdBergh-Hagen --> BH


##### MWSC+DIAS21+CG20+BICA19 matching

Databases from which the "OLD" databases were generated:

CG20 contains data from:
LIUPANG19 (35; LP_XXX), CASTRO-18-19-20 (550; UBC_XXX_), SIM19 (138, UPK_XXX)

DIAS21 contains data from:
LIUPANG19 (177; LP_XXX; DIAS21 lists clusters that are marked in LIUPAG19 as
"Class 3" or "need confirmation" but does not list some confirmed clusters in
LIUPAG19 such as LP_58), CASTRO-19-20 (287; UBC_XXX), SIM19 (127, UPK_XXX)

We remove the clusters that are present in the original databases from all
four DBs before processing. These are: LIUPANG19, CASTRO-18-19-20-22 (UBC),
and SIM19 (UPK) clusters.

There are 18 UBC clusters in BICA19 with 9 of them apparently associated to
a previously identified cluster: Alessi 161|UBC 3, TRSG 4|UBC 1,
Teutsch 179|UBC 9, AT 4|Alessi-Teutsch 4|Alessi 26|TRSG 3|UBC 12,
Alessi 190|UBC 5, Alessi 145|UBC 13, Renou 23|UBC 8, Alessi 94|UBC 2,
Teutsch 182|UBC 6. We remove the remaining 9.

After cross-matching these four databases (removing the UBC clusters) using
we have 4100 unique clusters.

```
DIAS21 1152
CG20 1294
MWSC 2858
BICA19_all_names 2895

Clusters in DIAS21 matched to clusters in CG20: 1093
Clusters in DIAS21 not matched to clusters in CG20: 59
Clusters in CG20 not matched to clusters in DIAS21: 201
Clusters in cross-matched catalog so far: 1353

Clusters in DIAS21+CG20 matched to clusters in MWSC: 1119
Clusters in DIAS21+CG20 not matched to clusters in MWSC: 234
Clusters in MWSC not matched to clusters in DIAS21+CG20: 1739
Clusters in cross-matched catalog so far: 3092

Clusters in DIAS21+CG20+MWSC matched to clusters in BICA19: 1904
Clusters in DIAS21+CG20+MWSC not matched to clusters in BICA19: 1188
Clusters in BICA19 not matched to clusters in DIAS21+CG20+MWSC: 1008
Clusters in final cross-matched catalog: 4100
```

For the Ryu clusters:
*The average angular radius of the clusters is 1'.31±0 60.
More specifically, 902 (98%) clusters are smaller than 3′, and
823 (89%) clusters are even smaller than 2′.* [Ryu & Lee (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...856..152R/abstract)

**UPDATE**
There are 1906 clusters with parallax data only from HAO21, which is of low
quality.

There are 1457 clusters with no proper motions data in any DB, which means that
there are 4137 clusters with PMs data in at least one of the DBs.
There are 1906 clusters with PMs only in HAO21, which is of low
quality.
This means that there are **3363** (1457+1906) clusters with either no or low
quality PMs data (ie: ~60% of the combined clusters).

There are 1457 clusters with neither parallax nor proper motions data: 707 are
Ryu, 388 are FSR, 51 are VVV, 46 are Teutsch, 39 are ESO, 30 are LS, 22 are
Bica, the remaining 174 are mixed. These same clusters contain no parallax and
no proper motions data separately (i.e.: there are no more clusters that do not
contain just parallax data, and no more clusters that do not contain just
proper motions data).


### Cluster definition

Consider a point c in an n-dimensional space and an integer value m>0. We define
a 'cluster' as the collection of m elements with the smallest n-dimensional
euclidean distance to c.




### IAU naming convention

https://cds.unistra.fr/Dic/iau-spec.html#S3.5.1




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
- []()

