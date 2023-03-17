
# NAME TO BE DETERMINED

This project began on June 2022.


## TODO

1. Clean up the `sc_full_cat` code
 1. ~~the OLD method works better but the NEW method correctly stores to which DB each cluster belongs to~~
 2. ~~Implement name matching prior to coordinates matching~~
 3. ~~finish re-formatting IO.readData() module~~
 4. ~~convert all input DBs to csv format~~
 5. Add databases listed below






### Databases

* [New Open-cluster Candidates Found in the Galactic Disk Using Gaia DR2/EDR3 Data, 2022](https://ui.adsabs.harvard.edu/abs/2022ApJS..260....8H/abstract): 541
* [Hunting for open clusters in Gaia EDR3: 628 new open clusters found with OCfinder, 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.118C/abstract): 628
* [Newly detected open clusters in the Galactic disk using Gaia EDR3, 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...660A...4H/abstract): 704
* [Updated parameters of 1743 open clusters based on Gaia DR2, 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.504..356D): 1743
* [Painting a portrait of the Galactic disc with its stellar clusters, 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...640A...1C) : 2017
* [Hunting for open clusters in Gaia DR2: 582 new open clusters .. Galactic disc, 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A..45C/abstract): 582 (only 570 in CDS table?)
* [A Catalog of Newly Identified Star Clusters in Gaia DR2, 2019](https://ui.adsabs.harvard.edu/abs/2019ApJS..245...32L/abstract): 76
* [A Multi-band Catalog of 10978 Star Clusters ... in the Milky Way, 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....157...12B/abstract): 3654

* [WEBDA](http://www.univie.ac.at/webda/): 1756
* [MWSC](https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html): 2858
* [OPENCLUST](https://heasarc.gsfc.nasa.gov/W3Browse/all/openclust.html): 2167



- 2022

* [A comprehensive photometric and kinematical characteristic of the newly discovered QCs clusters with Gaia EDR3](https://ui.adsabs.harvard.edu/abs/2022JApA...43...26E/abstract): 4

- 2021

* [Membership Lists for 431 Open Clusters in Gaia DR2 Using Extreme Deconvolution Gaussian Mixture Models](https://ui.adsabs.harvard.edu/abs/2021ApJ...923..129J/abstract): 11
* [New open clusters found by manual mining of data based on Gaia DR2](https://ui.adsabs.harvard.edu/abs/2021RAA....21..117C/abstract): 20
* [A catalogue of 74 new open clusters found in Gaia Data-Release 2](https://ui.adsabs.harvard.edu/abs/2021RAA....21...93H/abstract): 74
* [Discovery of Four New Clusters in the Cygnus Cloud](https://ui.adsabs.harvard.edu/abs/2021RAA....21...45Q/abstract): 4
* [New star clusters discovered towards the Galactic bulge direction using Gaia DR2](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502L..90F/abstract): 34

- 2020

* [Discovery and astrophysical properties of Galactic open clusters in dense stellar fields using Gaia DR2](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.2021F/abstract): 25
* [Sixteen Open Clusters Discovered with Sample-based Clustering Search of Gaia DR2](https://ui.adsabs.harvard.edu/abs/2020PASP..132c4502H/abstract): 16

- 2019

* [Gaia DR2 unravels incompleteness of nearby cluster population: new open clusters in the direction of Perseus](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.126C/abstract): 41
* [Three new Galactic star clusters discovered in the field of the open cluster NGC 5999 with Gaia DR2](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.5508F/abstract): 3

- pre-2018

* [Nine new open clusters within 500 pc from the Sun, 2016](https://ui.adsabs.harvard.edu/abs/2016A%26A...595A..22R/abstract): 9
* [Automated search for Galactic star clusters in large multiband surveys. I. Discovery of 15 new open clusters in the Galactic anticenter region, 2008](https://ui.adsabs.harvard.edu/abs/2008A%26A...486..771K/abstract): 15


### Cantat-Gaudin & Anders (2020)

1481 clusters (435833 stars with P>0.01); G_max=18

Average number of members

For members identified as those with P>0.5 (table1.csv):
 N<50 : 26.5%
 N<100: 57.1%
 N<200: 80.9%

For all the stars selected as members P>0.01 (members.csv):
 N<50 :  7.7% (115)
 N<100: 25.7% (380)
 N<200: 53.1% (787)
 N<300: 69.6% (1031)

Mean   = 294 stars per cluster (435833/1481)
Median = 184 stars per cluster
Max    = 3646 stars per cluster (NGC_7789)
Min    = 14 stars per cluster (DBSB_21)


### Cantat-Gaudin et al. (2020)

2017 clusters (234129 stars with P>0.7); parameters for 1867.
Mean = 125 stars per cluster (234129/1867)

1. What is the G limit used?
  The hard limit is 18, but there are a 15 stars up to ~19.6. These belong to
  Hyades/Melotte 25 (9) and Melotte 111 (6)
2. How were very large clusters processed?
  They were not processed with UPMASK. This is explained in Sect 2.1.

Large frames:

* r>30'  : 138 clusters
* d<1kpc : 245 clusters (there are 376 if we count those in the group above)

Sol processes 1634 clusters with r<30' & d>1kpc, I process 383 clusters: 138
with r>30' and 245 with d<1kpc.

Out of these 383 clusters there are 78 with one of these 2 flags:

1: r50>2
2: max_ra_dec>5

where 'r50' is the radius for half the members (according to CG2020)
'max_ra_dec'is the maximum distance between two members of this cluster in
either 'ra' or 'dec'.

```
Name          r50     max_ra_dec flag
-------------------------------------
Melotte_25    5.416   33.63      12
Melotte_20    2.027   25.45      12
Platais_8     2.716   23.96      12
Melotte_111   3.112   21.21      12
IC_2602       1.449   25.67      2
UPK_612       1.941   20.49      2
Alessi_3      1.541   16         2
Alessi_13     1.361   14.87      2
Platais_3     1.723   14.14      2
Alessi_9      1.73    13.86      2
Stock_2       1.03    13.67      2
NGC_6475      0.834   10.26      2
UPK_579       1.083   10.13      2
IC_2391       0.814   10.12      2
UPK_545       0.96    10.09      2
Ruprecht_98   0.835   9.96       2
BH_164        0.639   9.6        2
Collinder_463 0.608   9.53       2
UPK_640       1.213   9.41       2
UBC_11        1.034   9.18       2
Alessi_24     0.466   9.08       2
UPK_567       0.577   8.73       2
Alessi_5      0.542   8.62       2
Platais_9     1.595   8.59       2
Collinder_359 1.562   8.54       2
UPK_594       0.591   7.96       2
UPK_167       1.197   7.85       2
Melotte_22    1.274   7.81       2
UPK_350       0.838   7.62       2
Mamajek_4     1.404   7.44       2
NGC_2451A     1.146   7.39       2
UPK_560       1.005   7.21       2
UPK_166       0.767   7.07       2
BH_99         0.537   7.07       2
Collinder_135 1.138   6.95       2
ASCC_123      1.294   6.94       2
NGC_6991      0.555   6.84       2
Collinder_350 0.786   6.82       2
RSG_7         0.826   6.68       2
Stock_1       0.742   6.64       2
ASCC_73       0.988   6.6        2
ASCC_79       0.98    6.49       2
UPK_495       0.939   6.48       2
NGC_7092      1.003   6.36       2
NGC_3532      0.536   6.29       2
UPK_524       0.957   6.27       2
Trumpler_10   0.89    6.25       2
UBC_10a       0.512   6.24       2
UPK_303       0.734   6.19       2
RSG_8         1.015   6.13       2
UPK_168       0.781   6.12       2
LP_5          0.627   6.04       2
Gulliver_9    0.969   6.04       2
NGC_752       0.485   5.98       2
COIN-Gaia_13  1.003   5.93       2
UPK_533       1.016   5.87       2
Stock_23      0.754   5.87       2
UPK_305       1.597   5.86       2
UPK_535       1.136   5.63       2
Roslund_6     1.004   5.57       2
Ruprecht_147  0.652   5.56       2
NGC_1039      0.403   5.52       2
UPK_230       0.38    5.5        2
UPK_585       0.217   5.46       2
UPK_442       1.125   5.42       2
ASCC_127      0.627   5.41       2
NGC_3228      0.512   5.39       2
UPK_624       0.479   5.38       2
UPK_552       0.639   5.35       2
UBC_31        1.612   5.29       2
Alessi_2      0.554   5.27       2
NGC_6124      0.401   5.22       2
Collinder_132 1.333   5.21       2
UPK_296       0.808   5.19       2
UPK_526       0.809   5.14       2
ASCC_41       0.679   5.11       2
UPK_599       0.711   5.08       2
UBC_8         0.426   5.02       2
```

`Blanco_1` has a maximum longitude distance between members of 26.8, The '2'
flag is not raised because it wraps around 360 and 'max_ra_dec>100'. This
clusters is downloaded from vizier using Galactic center and a box size of 9º
and later on trimmed to a proper squared frame (I could not process it
properly with my own script to generate Gaia frames)

Two clusters in the subsample of 383 do not have a distance assigned:

```
Name,ra,dec,box,plx_d
---------------------
COIN-Gaia_37,45.377,58.329,4.225,nan
UPK_606,216.12,-46.39,3.23,nan
```

There are two clusters (Mamajek_1 & NGC_188) that are not in the subsample
of 383 clusters but have the '2' flag raised.


```
CG2020 2017 clusters
--------------------
  |  |
  |  |__ 1634 (smaller frames + larger distances)
  |     |
  |     |____ 426 (r50 >= 0.1)
  |     |    |
  |     |    |____ 26 (no distance assigned)
  |     |    |
  |     |    |____ 242 (d_pc <= 2 Kpc)
  |     |    |
  |     |    |____ 158 (d_pc > 2 Kpc)
  |     |
  |     |____ 623 (0.05 <= r50 < 0.1)
  |     |    |
  |     |    |____ 32 (no distance assigned)
  |     |    |
  |     |    |____ 174 (d_pc <= 2 Kpc)
  |     |    |
  |     |    |____ 417 (d_pc > 2 Kpc)
  |     |
  |     |____ 585 (r50 < 0.05)
  |          |
  |          |____ 90 (no distance assigned)
  |          |
  |          |____ 35 (d_pc <= 2 Kpc)
  |          |
  |          |____ 460 (d_pc > 2 Kpc)
  |
  |
  |___ 383 (largest frames + smaller distances)
      |
      |____ 2 (no distance assigned)
      |
      |____ 78 (large frame flags)
      |    |
      |    |____ 6 (size > 50Mb)
      |    |     MiniBatchKmeans? Voronoi?
      |    |
      |    |____ 32 (size > 10Mb)
      |    |     MiniBatchKmeans?
      |    |
      |    |____ 20 (size > 5Mb)
      |    |     Kmeans?
      |    |
      |    |____ 20 (size < 5Mb)
      |          GMM?
      |
      |____ 303 (no flags)
           |
           |____ 47 (size > 5Mb)
           |    |
           |    |____ 1 (size > 50Mb ; Gulliver_29)
           |    |     MiniBatchKmeans?
           |    |
           |    |____ 13 (size > 10Mb)
           |    |     MiniBatchKmeans?
           |    |
           |    |____ 33 (5Mb < size < 10Mb)
           |          Kmeans?
           |     
           |____ 256 (size < 5Mb)
                 GMM (+KDE,+GUMM) (~15 min each)
```

### Group of 256 large w no flags

This set was processed removing those stars below a given parallax threshold.
This is to improve the performance of pyUPMASK by removing the clearly non
members. The parallax limit was taken from the mean value assigned by CG2020,
and a region below that value was included (to contain the dispersed members)

The first run was processed using GMM. Visual inspection revealed that there
were XXX clusters that needed re-processing for various reasons.

**Notes on some clusters**

- UPK_433: shares central coordinates with Liu+Pang-2349



### Group of 47 large w no flags

xxx

#### Gulliver_29 (>50 Mb)

This is a very large frame (+70 Mb) from the no flag subset, for which I tested
several approaches.

1. MiniBatchKmeans (-KDE,-GUMM, 25 OL, 25 IL): reasonable res (~1 hr)
2. Voronoi (-KDE, -GUMM): bad res (~10 min)
3. MiniBatchKmeans (+KDE, +GUMM 50 internal runs, 5 OL, 10 IL): horrible res (~1 hr)
4. MiniBatchKmeans (-KDE, +GUMM 100 internal runs, 5 OL, 25 IL): horrible res (~35 min)
5. MiniBatchKmeans (+KDE, -GUMM, 5 OL, 50 IL): reasonable res, cleaner than 1 (~20 min)
6. MiniBatchKmeans (-KDE, -GUMM, 10 OL, 100 IL): reasonable res (~25 min)
7. Voronoi (+KDE, +GUMM, 100 IL): never finished (>1 hr)
8. Kmeans (-KDE, -GUMM, 10 OL, 50 IL): never finished (>1 hr)
8. KNN (+KDE, -GUMM, 50 IL, HDD not cdist): horrible res (~25 min)



#### Group of 13 large w no flags (10Mb < size < 50Mb)

Processed with: 3 dim, 25 OL, +resample, -PCA, -KDE, -GUMM, 50 IL, Kmeans.
It took ~75 min per cluster.

**Notes on some clusters**


#### Group of 33 large w no flags (5Mb < size < 10Mb)

Processed with: 3 dim, 25 OL, +resample, -PCA, -KDE, -GUMM, 50 IL, Kmeans.
It took ~45 min per cluster.


**Notes on some clusters**

- ASCC_32: shares field with Gulliver_21. Contains several overdensities in the
  VPD
- UPK_470: shares field with DBSB_3 (bad process in CG2020) and apparently
  another cluster (?) centered in the VPD at (-3.6, 4.2)
- UBC_26: shares field with Alessi_62



### Group of 78 (large frame flags)

**Notes on some clusters**

- RSG_7: shares a field 




### New fast method for membership estimation

Norm   : Normalized spaces in VPD and Plx
weight : weighting the distances by the uncertainties in PMs and parallax
GUMM   : Gaussian+Uniform mixture model used to filter likely non-members
std    : `C_thresh=np.median(C_S_field) + std * np.std(C_S_field)`


Runs 1 through 8 all have their VPD+Plx dimensions normalized.
Results:
- no GUMM means a larger dispersion in coordinates (as expected)
- weighting by errors increases the percentage lost, mostly in the +16 range
- std=2 makes very little difference

Runs 9 through 16 don't normalize the VPD+Plx dimensions.
Results:
- no GUMM means a larger dispersion in coordinates (as expected)
- larger dispersion in Plx compared to runs 1-8
- more CG members are on average recovered compared to runs 1-8, mainly in
  the +16 range
- Runs 15 and 13 miss 100% of CG2020 members for some clusters

Runs 17 through 24 are processed with a smaller `prob_cut` for the GUMM.
Results:
- No visible improvement

Runs 25 through 32 are processed using the averaged GUMM probabilities (ie: no 
probability cut)
Results:
- Bad results, increases the missed percentage in both ranges (-16 and +16)

Runs 33 through 40 skipped.

Runs 41 through 43 apply the second new method where coordinates distances are
also used to assign membership likelihoods. Runs 41 and 42 gave poor results
so I removed them. Run 43 is useful but I have no way to actually select the
most probable cluster members.

Runs 44 to 49 apply different outlier removal methods: isol=IsolationForest,
LOF=LocalOutlierFactor, NSTD=sigma clipping. The sigma clipping is the fastest
and the one that gives best results.

Runs 50 to 57 refine on run 49 which was the best of the previous batch.
Results:
- Runs with std=2 miss a lot of CG2020 members, run 51 is the best

Runs 58 to 65 are equivalent to 50-57 but they use the VPD+Plx space to reject
outliers instead of the xy space.
Results:
- Runs with std=2 miss a lot of CG2020 members,

Runs 66-67 use larger clusters in the Rk block, `2*N_membs` and `3*N_membs`
respectively.
Results:
- The probability cut value becomes erratic

Runs 68 through 71 don't apply any kind of GUMM, instead the classification
probability is obtained and finally averaged. Tested this with the version of
code used for the runs 51 and later.
Results:
- The KDE classifier makes the process much slower. Using all 5 dimensions it
introduces noise and pushes probabilities to more extreme values with no
visible reward
- Tested the LogisticRegression but got bad results
- GaussianNB+isotonic losses members for some clusters
- GaussianNB+sigmoid better at recovering members but includes many non-members.
  Also the probabilities tend to be lowered.

Run 72 uses ASteCA's Bayesian method to assign probabilities.
Results:
- The probabilities are lowered and no real improvement is gained.

Run 73 normalizes the VPD+Plx data before estimating the distances to the
centers.
Results:
- Lower probabilities compared to GMM
- More noise in the CMD
- Low mag stars can be lost
- Dispersion in the coordinates space

Run 74 adds to 73 the estimation of the VPD+Plx center outside of the
resampling for block instead of inside
Results:
- Lower probabilities compared to GMM
- Low mag stars can be lost
- Dispersion in the coordinates space



These use `MinMaxScaler` in xy and `C_thresh` break + `N_break` reset
Run 75: normalize, C_thresh 1.42
Run 76: non normalize, C_thresh 1.68  <-- Best
Run 77: normalize, C_thresh 1.68
Run 78: non normalize, C_thresh 1.42

Results (GMM cut P=0.7, fastMP cut P=0.5):
- 75: low magnitude GMM members are not recovered
- 76: recovers low magnitude members, adds noise to CMD
- 77: some low magnitude GMM members are lost
- 78: some low magnitude GMM members are lost (juchert_18), noise in CMD
(ngc_2259)

These use `MinMaxScaler` and `C_S_field` break
Run 79: normalize
Run 80: non normalize  <-- Best

So far 76 is better than 80 because it adds less noise


These use no `MinMaxScaler`, outlier reject, and `C_thresh` break with no
`N_break` reset
Run 83: normalize, C_thresh 1.42
Run 84: non normalize, C_thresh 1.68  <-- Best
Run 85: normalize, C_thresh 1.68
Run 86: non normalize, C_thresh 1.42

So far 84 is better than 76 because it adds less noise
84: juchert_18 failed

No `MinMaxScaler`, outlier reject, and `C_thresh` break with `N_break` reset
Run 87: non normalize, C_thresh 1.68

Runs 76, 84, and 87 are very similar with the exception of juschert_18 where
only 76 gives reasonable results.

No `MinMaxScaler`, outlier reject, and `C_thresh` break with `N_break` reset
Run 88: non normalize, `C_thresh` 1.68, add extra `C_thresh` reject at P>0.5
Results are bad, many members are lost, even bright ones