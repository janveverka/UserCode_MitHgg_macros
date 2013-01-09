#ifndef MODEL_CONSTANTS
#define MODEL_CONSTANTS

#include "TString.h"

void initMassArrayMap(int numMasses, double* masses, std::map<double,int>& map) {
  
  for(int iMass=0; iMass<numMasses; ++iMass)
    map.insert(std::pair<double,int>(masses[iMass],iMass));

  return;
}

//---making cross section and branching ratio as funciton of mass---
const int numsmpoints = 81;
//mass 
double smmasses[numsmpoints] = {110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,140.5,141.0,141.5,142.0,142.5,143.0,143.5,144.0,144.5,145.0,145.5,146.0,146.5,147.0,147.5,148.0,148.5,149.0,149.5,150.0};  


//branching ratio
double smbr[numsmpoints] = {0.00197,0.00199,0.002,0.00202,0.00204,0.00205,0.00207,0.00209,0.0021,0.00212,0.00213,0.00215,0.00216,0.00217,0.00218,0.0022,0.00221,0.00222,0.00223,0.00224,0.00225,0.00226,0.00226,0.00227,0.00228,0.00228,0.00228,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00229,0.00228,0.00228,0.00227,0.00227,0.00226,0.00225,0.00224,0.00223,0.00222,0.00221,0.00219,0.00218,0.00217,0.00215,0.00213,0.00212,0.0021,0.00208,0.00206,0.00204,0.00202,0.002,0.00198,0.00196,0.00193,0.00191,0.00189,0.001865,0.00184,0.00181,0.00178,0.001755,0.00173,0.0017,0.00167,0.001645,0.00162,0.00159,0.00156,0.001525,0.00149,0.00146,0.00143,0.001395,0.00136};
double ffbr[numsmpoints] = {0.059540,0.056510,0.053660,0.050980,0.048460,0.046090,0.043860,0.041760,0.039780,0.037910,0.036150,0.034490,0.032920,0.031430,0.030030,0.028710,0.027460,0.026270,0.025150,0.024080,0.023070,0.022120,0.021210,0.020340,0.019520,0.018740,0.018000,0.017300,0.016630,0.015990,0.015380,0.014800,0.014240,0.013710,0.013210,0.012720,0.012260,0.011820,0.011400,0.010990,0.010610,0.010240,0.009880,0.009539,0.009212,0.008897,0.008595,0.008305,0.008026,0.007758,0.007500,0.007251,0.007012,0.006781,0.006558,0.006344,0.006137,0.005937,0.005744,0.005557,0.005377,0.005202,0.005034,0.004870,0.004711,0.004558,0.004409,0.004264,0.004124,0.003987,0.003855,0.003725,0.003600,0.003477,0.003358,0.003241,0.003128,0.003016,0.002908,0.002801,0.002697};
double sm4br[numsmpoints] = {4.40E-005,0.000043685,0.00004337,0.000043055,0.00004274,0.000042425,0.00004211,0.000041795,0.00004148,0.000041165,0.00004085,0.000040535,0.00004022,0.000039905,0.00003959,0.000039275,0.00003896,0.000038645,0.00003833,0.000038015,3.77E-005,0.00003717,0.00003664,0.00003611,0.00003558,0.00003505,0.00003452,0.00003399,0.00003346,0.00003293,0.0000324,0.00003187,0.00003134,0.00003081,0.00003028,0.00002975,0.00002922,0.00002869,0.00002816,0.00002763,2.71E-005,0.000026395,0.00002569,0.000024985,0.00002428,0.000023575,0.00002287,0.000022165,0.00002146,0.000020755,0.00002005,0.000019345,0.00001864,0.000017935,0.00001723,0.000016525,0.00001582,0.000015115,0.00001441,0.000013705,1.30E-005,0.000012421,0.000011842,0.000011263,0.000010684,0.000010105,0.000009526,0.000008947,0.000008368,0.000007789,0.00000721,0.000006631,0.000006052,0.000005473,0.000004894,0.000004315,0.000003736,0.000003157,0.000002578,0.000001999,1.42E-006};

//cross section for different production mechanism (OLD NUMBERS)
double gghxsec[numsmpoints] = {
  25.0119,
  24.7918,
  24.5743,
  24.3597,
  24.1486,
  23.9401,
  23.7343,
  23.5312,
  23.3306,
  23.1325,
  22.9369,
  22.7442,
  22.5534,
  22.3651,
  22.1790,
  21.9949,
  21.8134,
  21.6341,
  21.4570,
  21.2820,
  21.1090,
  20.9381,
  20.7692,
  20.6023,
  20.4373,
  20.2742,
  20.1131,
  19.9538,
  19.7963,
  19.6407,
  19.4868,
  19.3347,
  19.1843,
  19.0357,
  18.8887,
  18.7433,
  18.5996,
  18.4574,
  18.3169,
  18.1778,
  18.0403,
  17.9043,
  17.7698,
  17.6368,
  17.5053,
  17.3752,
  17.2464,
  17.1191,
  16.9930,
  16.8683,
  16.7448,
  16.6226,
  16.5017,
  16.3820,
  16.2635,
  16.1463,
  16.0303,
  15.9155,
  15.8019,
  15.6894,
  15.5782,
  (15.5782+15.3592)/2.,
  15.3592,
  (15.3592+15.1446)/2.,
  15.1446,
  (15.1446+14.9342)/2.,
  14.9342,
  (14.9342+14.7278)/2.,
  14.7278,
  (14.7278+14.5253)/2.,
  14.5253,
  (14.5253+14.3267)/2.,
  14.3267,
  (14.3267+14.1318)/2.,
  14.1318,
  (14.1318+13.9404)/2.,
  13.9404,
  (13.9404+13.7522)/2.,
  13.7522,
  (13.7522+13.5672)/2.,
  13.5672
};

double vbfxsec[numsmpoints] = {
  1.791,
  1.783,
  1.775,
  1.766,
  1.758,
  1.750,
  1.742,
  1.733,
  1.725,
  1.717,
  1.709,
  1.701,
  1.693,
  1.686,
  1.678,
  1.670,
  1.661,
  1.654,
  1.647,
  1.639,
  1.632,
  1.624,
  1.617,
  1.609,
  1.602,
  1.595,
  1.588,
  1.580,
  1.573,
  1.566,
  1.559,
  1.552,
  1.544,
  1.539,
  1.531,
  1.524,
  1.517,
  1.511,
  1.504,
  1.497,
  1.490,
  1.483,
  1.477,
  1.470,
  1.463,
  1.458,
  1.451,
  1.444,
  1.439,
  1.432,
  1.425,
  1.419,
  1.413,
  1.407,
  1.401,
  1.395,
  1.388,
  1.382,
  1.376,
  1.370,
  1.365,
  (1.365+1.352)/2.,
  1.352,
  (1.352+1.341)/2.,
  1.341,
  (1.341+1.329)/2.,
  1.329,
  (1.329+1.317)/2.,
  1.317,
  (1.317+1.306)/2.,
  1.306,
  (1.306+1.295)/2.,
  1.295,
  (1.295+1.284)/2.,
  1.284,
  (1.284+1.272)/2.,
  1.272,
  (1.272+1.261)/2.,
  1.261,
  (1.261+1.251)/2.,
  1.251
};

double whxsec[numsmpoints] = {
  1.060,
  1.045,
  1.030,
  1.015,
  0.9998,
  0.9852,
  0.9709,
  0.9570,
  0.9432,
  0.9297,
  0.9165,
  0.9035,
  0.8907,
  0.8782,
  0.8659,
  0.8538,
  0.8420,
  0.8303,
  0.8187,
  0.8075,
  0.7966,
  0.7859,
  0.7753,
  0.7649,
  0.7547,
  0.7446,
  0.7347,
  0.7249,
  0.7154,
  0.7060,
  0.6966,
  0.6873,
  0.6782,
  0.6691,
  0.6602,
  0.6515,
  0.6429,
  0.6344,
  0.6260,
  0.6177,
  0.6095,
  0.6015,
  0.5936,
  0.5859,
  0.5783,
  0.5708,
  0.5634,
  0.5562,
  0.5491,
  0.5420,
  0.5351,
  0.5283,
  0.5215,
  0.5149,
  0.5084,
  0.5020,
  0.4956,
  0.4894,
  0.4833,
  0.4772,
  0.4713,
  (0.4713 + 0.4597)/2. ,
  0.4597,
  (0.4597 + 0.4484)/2.,
  0.4484,
  (0.4484 + 0.4375)/2.,
  0.4375,
  (0.4375 + 0.4268)/2.,
  0.4268,
  (0.4268 + 0.4164)/2. ,
  0.4164,
  (0.4164 + 0.4062)/2.,
  0.4062,
  (0.4062 + 0.3963)/2.,
  0.3963,
  (0.3963 + 0.3867)/2.,
  0.3867,
  (0.3867 + 0.3773)/2.,
  0.3773,
  (0.3773 +0.3681)/2.,
  0.3681
};
  
//double whxsec[numsmpoints] = {0.8754,0.8623,0.8495,0.8368,0.8244,0.8122,0.8003,0.7885,0.7770,0.7657,0.7546,0.7439,0.7333,0.7230,0.7129,0.7030,0.6933,0.6837,0.6744,0.6651,0.6561,0.6472,0.6384,0.6297,0.6212,0.6129,0.6046,0.5965,0.5885,0.5806,0.5729,0.5652,0.5576,0.5501,0.5428,0.5355,0.5284,0.5213,0.5144,0.5075,0.5008,0.4942,0.4877,0.4813,0.4749,0.4687,0.4626,0.4566,0.4506,0.4448,0.4390,0.4333,0.4277,0.4221,0.4167,0.4113,0.4060,0.4008,0.3957,0.3907,0.3857,0.3809,0.3761,0.3715,0.3669,0.3624,0.3579,0.3535,0.3491,0.3449,0.3406,0.3364,0.3321,0.3280,0.3238,0.3198,0.3157,0.3118,0.3078,0.3040,0.3001};

double zhxsec[numsmpoints] = {
  0.5869,
  0.5788,
  0.5708,
  0.5629,
  0.5552,
  0.5476,
  0.5402,
  0.5329,
  0.5258,
  0.5187,
  0.5117,
  0.5049,
  0.4981,
  0.4916,
  0.4850,
  0.4787,
  0.4724,
  0.4662,
  0.4602,
  0.4542,
  0.4483,
  0.4426,
  0.4368,
  0.4312,
  0.4257,
  0.4203,
  0.4150,
  0.4096,
  0.4044,
  0.3993,
  0.3943,
  0.3893,
  0.3843,
  0.3794,
 	0.3746,
 	0.3699,
 	0.3652,
 	0.3606,
 	0.3561,
 	0.3516,
 	0.3473,
 	0.3430,
 	0.3388,
 	0.3347,
 	0.3306,
 	0.3266,
 	0.3226,
 	0.3188,
 	0.3149,
 	0.3112,
 	0.3074,
 	0.3038,
 	0.3001,
 	0.2966,
 	0.2930,
 	0.2895,
 	0.2861,
 	0.2827,
 	0.2793,
 	0.2760,
 	0.2728,
(0.2728+0.2664)/2.,
 	0.2664,
(0.2664+0.2601)/2.,
 	0.2601,
(0.2601+0.2541)/2.,
 	0.2541,
(0.2541+0.2482)/2.,
 	0.2482,
(0.2482+0.2424)/2.,
 	0.2424,
(0.2424+0.2368)/2.,
 	0.2368,
(0.2368+0.2314)/2.,
 	0.2314,
(0.2314+0.2261)/2.,
 	0.2261,
(0.2261+0.2159)/2.,
 	0.2209,
(0.2209+0.2159)/2.,
 	0.2159
};

//double zhxsec[numsmpoints] = {0.4721,0.4655,0.4589,0.4525,0.4462,0.4400,0.4340,0.4280,0.4221,0.4164,0.4107,0.4052,0.3998,0.3945,0.3893,0.3842,0.3791,0.3742,0.3693,0.3645,0.3598,0.3551,0.3505,0.3459,0.3414,0.3370,0.3326,0.3283,0.3241,0.3199,0.3158,0.3117,0.3077,0.3038,0.2999,0.2961,0.2923,0.2886,0.2849,0.2813,0.2778,0.2743,0.2709,0.2675,0.2642,0.2609,0.2577,0.2545,0.2514,0.2483,0.2453,0.2423,0.2393,0.2364,0.2336,0.2307,0.2279,0.2252,0.2225,0.2198,0.2172,0.2147,0.2121,0.2096,0.2071,0.2047,0.2023,0.2000,0.1976,0.1953,0.1930,0.1907,0.1884,0.1862,0.1840,0.1818,0.1796,0.1775,0.1754,0.1734,0.1713};

double tthxsec[numsmpoints] = {
 	0.1887,
 	0.1863,
 	0.1840,
 	0.1816,
 	0.1793,
 	0.1771,
 	0.1748,
 	0.1727,
 	0.1705,
 	0.1684,
 	0.1663,
 	0.1642,
 	0.1622,
 	0.1602,
 	0.1582,
 	0.1563,
 	0.1543,
 	0.1524,
 	0.1506,
 	0.1487,
 	0.1470,
 	0.1452,
 	0.1434,
 	0.1417,
 	0.1400,
 	0.1383,
 	0.1366,
	0.1350,
 	0.1334,
 	0.1318,
 	0.1302,
 	0.1287,
	0.1271,
 	0.1256,
 	0.1241,
 	0.1227,
 	0.1212,
 	0.1198,
 	0.1184,
 	0.1170,
 	0.1157,
 	0.1143,
 	0.1130,
 	0.1117,
 	0.1104,
 	0.1091,
 	0.1079,
 	0.1067,
 	0.1054,
 	0.1042,
 	0.1031,
 	0.1019,
 	0.1007,
 	0.09961,
 	0.09849,
 	0.09738,
 	0.09629,
 	0.09521,
 	0.09415,
 	0.09310,
 	0.09207,
(0.09207+0.09004)/2.,
 	0.09004,
(0.09004+0.08807)/2.,
 	0.08807,
(0.08807+0.08615)/2.,
 	0.08615,
(0.08615+0.08427)/2.,
 	0.08427,
(0.08427+0.08246)/2.,
 	0.08246,
(0.08246+0.08068)/2.,
 	0.08068,
(0.08068+0.07895)/2.,
 	0.07895,
(0.07895+0.07727)/2.,
 	0.07727,
(0.07727+0.07563)/2.,
 	0.07563,
(0.07563+0.07403)/2.,
 	0.07403
};

//double tthxsec[numsmpoints] = {0.12570,0.12410,0.12250,0.12090,0.11940,0.11790,0.11640,0.11490,0.11340,0.11200,0.11060,0.10920,0.10780,0.10650,0.10510,0.10380,0.10250,0.10130,0.10000,0.09878,0.09756,0.09636,0.09518,0.09402,0.09287,0.09174,0.09063,0.08954,0.08846,0.08739,0.08634,0.08530,0.08428,0.08327,0.08227,0.08129,0.08032,0.07937,0.07842,0.07750,0.07658,0.07568,0.07479,0.07391,0.07304,0.07219,0.07135,0.07052,0.06970,0.06890,0.06810,0.06731,0.06654,0.06577,0.06502,0.06428,0.06355,0.06282,0.06211,0.06141,0.06072,0.06005,0.05937,0.05872,0.05807,0.05744,0.05680,0.05618,0.05556,0.05496,0.05435,0.05376,0.05316,0.05258,0.05200,0.05144,0.05087,0.05032,0.04976,0.04923,0.04869};

double sm4gghxsec[numsmpoints] = {199,197.3,195.6,193.9,192.2,190.5,188.8,187.1,185.4,183.7,182,180.3,178.6,176.9,175.2,173.5,171.8,170.1,168.4,166.7,165,163.65,162.3,160.95,159.6,158.25,156.9,155.55,154.2,152.85,151.5,150.15,148.8,147.45,146.1,144.75,143.4,142.05,140.7,139.35,138,136.95,135.9,134.85,133.8,132.75,131.7,130.65,129.6,128.55,127.5,126.45,125.4,124.35,123.3,122.25,121.2,120.15,119.1,118.05,117,116,115,114,113,112,111,110,109,108,107,106.22,105.44,104.66,103.88,103.1,102.32,101.54,100.76,99.98,99.2};

double wzhxsec[numsmpoints];

double nullarray[numsmpoints];

void initProcessCSArrayMap(std::map<TString,double*>& processCrossSectionMap) {
  
  // FIX-ME: this is a hack for now, move to use only maps ....
  std::map<double,double> XSectionMap_ggh;
  std::map<double,double> XSectionMap_vbf;
  std::map<double,double> XSectionMap_wh;
  std::map<double,double> XSectionMap_zh;
  std::map<double,double> XSectionMap_tth;
  
  XSectionMap_ggh.insert(std::pair<double,double>(110.0,25.04));
  XSectionMap_ggh.insert(std::pair<double,double>(110.5,24.82));
  XSectionMap_ggh.insert(std::pair<double,double>(111.0,24.60));
  XSectionMap_ggh.insert(std::pair<double,double>(111.5,24.39));
  XSectionMap_ggh.insert(std::pair<double,double>(112.0,24.18));
  XSectionMap_ggh.insert(std::pair<double,double>(112.5,23.97));
  XSectionMap_ggh.insert(std::pair<double,double>(113.0,23.76));
  XSectionMap_ggh.insert(std::pair<double,double>(113.5,23.56));
  XSectionMap_ggh.insert(std::pair<double,double>(114.0,23.36));
  XSectionMap_ggh.insert(std::pair<double,double>(114.5,23.16));
  XSectionMap_ggh.insert(std::pair<double,double>(115.0,22.96));
  XSectionMap_ggh.insert(std::pair<double,double>(115.5,22.77));
  XSectionMap_ggh.insert(std::pair<double,double>(116.0,22.58));
  XSectionMap_ggh.insert(std::pair<double,double>(116.5,22.39));
  XSectionMap_ggh.insert(std::pair<double,double>(117.0,22.20));
  XSectionMap_ggh.insert(std::pair<double,double>(117.5,22.02));
  XSectionMap_ggh.insert(std::pair<double,double>(118.0,21.84));
  XSectionMap_ggh.insert(std::pair<double,double>(118.5,21.66));
  XSectionMap_ggh.insert(std::pair<double,double>(119.0,21.48));
  XSectionMap_ggh.insert(std::pair<double,double>(119.5,21.31));
  XSectionMap_ggh.insert(std::pair<double,double>(120.0,21.13));
  XSectionMap_ggh.insert(std::pair<double,double>(120.5,20.96));
  XSectionMap_ggh.insert(std::pair<double,double>(121.0,20.80));
  XSectionMap_ggh.insert(std::pair<double,double>(121.5,20.63));
  XSectionMap_ggh.insert(std::pair<double,double>(122.0,20.48));
  XSectionMap_ggh.insert(std::pair<double,double>(122.5,20.31));
  XSectionMap_ggh.insert(std::pair<double,double>(123.0,20.15));
  XSectionMap_ggh.insert(std::pair<double,double>(123.5,19.99));
  XSectionMap_ggh.insert(std::pair<double,double>(124.0,19.83));
  XSectionMap_ggh.insert(std::pair<double,double>(124.5,19.68));
  XSectionMap_ggh.insert(std::pair<double,double>(125.0,19.52));
  XSectionMap_ggh.insert(std::pair<double,double>(125.5,19.37));
  XSectionMap_ggh.insert(std::pair<double,double>(126.0,19.22));
  XSectionMap_ggh.insert(std::pair<double,double>(126.5,19.07));
  XSectionMap_ggh.insert(std::pair<double,double>(127.0,18.92));
  XSectionMap_ggh.insert(std::pair<double,double>(127.5,18.78));
  XSectionMap_ggh.insert(std::pair<double,double>(128.0,18.63));
  XSectionMap_ggh.insert(std::pair<double,double>(128.5,18.49));
  XSectionMap_ggh.insert(std::pair<double,double>(129.0,18.35));
  XSectionMap_ggh.insert(std::pair<double,double>(129.5,18.21));
  XSectionMap_ggh.insert(std::pair<double,double>(130.0,18.07));
  XSectionMap_ggh.insert(std::pair<double,double>(130.5,17.94));
  XSectionMap_ggh.insert(std::pair<double,double>(131.0,17.81));
  XSectionMap_ggh.insert(std::pair<double,double>(131.5,17.68));
  XSectionMap_ggh.insert(std::pair<double,double>(132.0,17.55));
  XSectionMap_ggh.insert(std::pair<double,double>(132.5,17.42));
  XSectionMap_ggh.insert(std::pair<double,double>(133.0,17.29));
  XSectionMap_ggh.insert(std::pair<double,double>(133.5,17.16));
  XSectionMap_ggh.insert(std::pair<double,double>(134.0,17.04));
  XSectionMap_ggh.insert(std::pair<double,double>(134.5,16.92));
  XSectionMap_ggh.insert(std::pair<double,double>(135.0,16.79));
  XSectionMap_ggh.insert(std::pair<double,double>(135.5,16.67));
  XSectionMap_ggh.insert(std::pair<double,double>(136.0,16.55));
  XSectionMap_ggh.insert(std::pair<double,double>(136.5,16.43));
  XSectionMap_ggh.insert(std::pair<double,double>(137.0,16.31));
  XSectionMap_ggh.insert(std::pair<double,double>(137.5,16.20));
  XSectionMap_ggh.insert(std::pair<double,double>(138.0,16.08));
  XSectionMap_ggh.insert(std::pair<double,double>(138.5,15.96));
  XSectionMap_ggh.insert(std::pair<double,double>(139.0,15.85));
  XSectionMap_ggh.insert(std::pair<double,double>(139.5,15.74));
  XSectionMap_ggh.insert(std::pair<double,double>(140.0,15.63));
  XSectionMap_ggh.insert(std::pair<double,double>(141.0,15.42));
  XSectionMap_ggh.insert(std::pair<double,double>(142.0,15.20));
  XSectionMap_ggh.insert(std::pair<double,double>(143.0,15.00));
  XSectionMap_ggh.insert(std::pair<double,double>(144.0,14.79));
  XSectionMap_ggh.insert(std::pair<double,double>(145.0,14.59));
  XSectionMap_ggh.insert(std::pair<double,double>(146.0,14.40));
  XSectionMap_ggh.insert(std::pair<double,double>(147.0,14.21));
  XSectionMap_ggh.insert(std::pair<double,double>(148.0,14.02));
  XSectionMap_ggh.insert(std::pair<double,double>(149.0,13.83));
  XSectionMap_ggh.insert(std::pair<double,double>(150.0,13.65));
    
    
  XSectionMap_vbf.insert(std::pair<double,double>(110.0,1.809));
  XSectionMap_vbf.insert(std::pair<double,double>(110.5,1.799));
  XSectionMap_vbf.insert(std::pair<double,double>(111.0,1.791));
  XSectionMap_vbf.insert(std::pair<double,double>(111.5,1.784));
  XSectionMap_vbf.insert(std::pair<double,double>(112.0,1.780));
  XSectionMap_vbf.insert(std::pair<double,double>(112.5,1.771));
  XSectionMap_vbf.insert(std::pair<double,double>(113.0,1.764));
  XSectionMap_vbf.insert(std::pair<double,double>(113.5,1.753));
  XSectionMap_vbf.insert(std::pair<double,double>(114.0,1.743));
  XSectionMap_vbf.insert(std::pair<double,double>(114.5,1.735));
  XSectionMap_vbf.insert(std::pair<double,double>(115.0,1.729));
  XSectionMap_vbf.insert(std::pair<double,double>(115.5,1.719));
  XSectionMap_vbf.insert(std::pair<double,double>(116.0,1.714));
  XSectionMap_vbf.insert(std::pair<double,double>(116.5,1.704));
  XSectionMap_vbf.insert(std::pair<double,double>(117.0,1.699));
  XSectionMap_vbf.insert(std::pair<double,double>(117.5,1.688));
  XSectionMap_vbf.insert(std::pair<double,double>(118.0,1.683));
  XSectionMap_vbf.insert(std::pair<double,double>(118.5,1.675));
  XSectionMap_vbf.insert(std::pair<double,double>(119.0,1.666));
  XSectionMap_vbf.insert(std::pair<double,double>(119.5,1.659));
  XSectionMap_vbf.insert(std::pair<double,double>(120.0,1.649));
  XSectionMap_vbf.insert(std::pair<double,double>(120.5,1.643));
  XSectionMap_vbf.insert(std::pair<double,double>(121.0,1.636));
  XSectionMap_vbf.insert(std::pair<double,double>(121.5,1.631));
  XSectionMap_vbf.insert(std::pair<double,double>(122.0,1.623));
  XSectionMap_vbf.insert(std::pair<double,double>(122.5,1.615));
  XSectionMap_vbf.insert(std::pair<double,double>(123.0,1.608));
  XSectionMap_vbf.insert(std::pair<double,double>(123.5,1.598));
  XSectionMap_vbf.insert(std::pair<double,double>(124.0,1.595));
  XSectionMap_vbf.insert(std::pair<double,double>(124.5,1.587));
  XSectionMap_vbf.insert(std::pair<double,double>(125.0,1.578));
  XSectionMap_vbf.insert(std::pair<double,double>(125.5,1.573));
  XSectionMap_vbf.insert(std::pair<double,double>(126.0,1.568));
  XSectionMap_vbf.insert(std::pair<double,double>(126.5,1.558));
  XSectionMap_vbf.insert(std::pair<double,double>(127.0,1.552));
  XSectionMap_vbf.insert(std::pair<double,double>(127.5,1.524));
  XSectionMap_vbf.insert(std::pair<double,double>(128.0,1.540));
  XSectionMap_vbf.insert(std::pair<double,double>(128.5,1.531));
  XSectionMap_vbf.insert(std::pair<double,double>(129.0,1.525));
  XSectionMap_vbf.insert(std::pair<double,double>(129.5,1.513));
  XSectionMap_vbf.insert(std::pair<double,double>(130.0,1.511));
  XSectionMap_vbf.insert(std::pair<double,double>(130.5,1.504));
  XSectionMap_vbf.insert(std::pair<double,double>(131.0,1.497));
  XSectionMap_vbf.insert(std::pair<double,double>(131.5,1.492));
  XSectionMap_vbf.insert(std::pair<double,double>(132.0,1.485));
  XSectionMap_vbf.insert(std::pair<double,double>(132.5,1.479));
  XSectionMap_vbf.insert(std::pair<double,double>(133.0,1.473));
  XSectionMap_vbf.insert(std::pair<double,double>(133.5,1.466));
  XSectionMap_vbf.insert(std::pair<double,double>(134.0,1.462));
  XSectionMap_vbf.insert(std::pair<double,double>(134.5,1.455));
  XSectionMap_vbf.insert(std::pair<double,double>(135.0,1.448));
  XSectionMap_vbf.insert(std::pair<double,double>(135.5,1.444));
  XSectionMap_vbf.insert(std::pair<double,double>(136.0,1.436));
  XSectionMap_vbf.insert(std::pair<double,double>(136.5,1.429));
  XSectionMap_vbf.insert(std::pair<double,double>(137.0,1.423));
  XSectionMap_vbf.insert(std::pair<double,double>(137.5,1.417));
  XSectionMap_vbf.insert(std::pair<double,double>(138.0,1.412));
  XSectionMap_vbf.insert(std::pair<double,double>(138.5,1.407));
  XSectionMap_vbf.insert(std::pair<double,double>(139.0,1.400));
  XSectionMap_vbf.insert(std::pair<double,double>(139.5,1.396));
  XSectionMap_vbf.insert(std::pair<double,double>(140.0,1.389));
  XSectionMap_vbf.insert(std::pair<double,double>(141.0,1.377));
  XSectionMap_vbf.insert(std::pair<double,double>(142.0,1.365));
  XSectionMap_vbf.insert(std::pair<double,double>(143.0,1.354));
  XSectionMap_vbf.insert(std::pair<double,double>(144.0,1.344));
  XSectionMap_vbf.insert(std::pair<double,double>(145.0,1.333));
  XSectionMap_vbf.insert(std::pair<double,double>(146.0,1.321));
  XSectionMap_vbf.insert(std::pair<double,double>(147.0,1.311));
  XSectionMap_vbf.insert(std::pair<double,double>(148.0,1.302));
  XSectionMap_vbf.insert(std::pair<double,double>(149.0,1.291));
  XSectionMap_vbf.insert(std::pair<double,double>(150.0,1.280));
    
  XSectionMap_wh.insert(std::pair<double,double>(110.0,1.060));
  XSectionMap_wh.insert(std::pair<double,double>(110.5,1.045));
  XSectionMap_wh.insert(std::pair<double,double>(111.0,1.030));
  XSectionMap_wh.insert(std::pair<double,double>(111.5,1.015));
  XSectionMap_wh.insert(std::pair<double,double>(112.0,0.9998));
  XSectionMap_wh.insert(std::pair<double,double>(112.5,0.9852));
  XSectionMap_wh.insert(std::pair<double,double>(113.0,0.9709));
  XSectionMap_wh.insert(std::pair<double,double>(113.5,0.9570));
  XSectionMap_wh.insert(std::pair<double,double>(114.0,0.9432));
  XSectionMap_wh.insert(std::pair<double,double>(114.5,0.9297));
  XSectionMap_wh.insert(std::pair<double,double>(115.0,0.9165));
  XSectionMap_wh.insert(std::pair<double,double>(115.5,0.9035));
  XSectionMap_wh.insert(std::pair<double,double>(116.0,0.8907));
  XSectionMap_wh.insert(std::pair<double,double>(116.5,0.8782));
  XSectionMap_wh.insert(std::pair<double,double>(117.0,0.8659));
  XSectionMap_wh.insert(std::pair<double,double>(117.5,0.8538));
  XSectionMap_wh.insert(std::pair<double,double>(118.0,0.8420));
  XSectionMap_wh.insert(std::pair<double,double>(118.5,0.8303));
  XSectionMap_wh.insert(std::pair<double,double>(119.0,0.8187));
  XSectionMap_wh.insert(std::pair<double,double>(119.5,0.8075));
  XSectionMap_wh.insert(std::pair<double,double>(120.0,0.7966));
  XSectionMap_wh.insert(std::pair<double,double>(120.5,0.7859));
  XSectionMap_wh.insert(std::pair<double,double>(121.0,0.7753));
  XSectionMap_wh.insert(std::pair<double,double>(121.5,0.7649));
  XSectionMap_wh.insert(std::pair<double,double>(122.0,0.7547));
  XSectionMap_wh.insert(std::pair<double,double>(122.5,0.7446));
  XSectionMap_wh.insert(std::pair<double,double>(123.0,0.7347));
  XSectionMap_wh.insert(std::pair<double,double>(123.5,0.7249));
  XSectionMap_wh.insert(std::pair<double,double>(124.0,0.7154));
  XSectionMap_wh.insert(std::pair<double,double>(124.5,0.7060));
  XSectionMap_wh.insert(std::pair<double,double>(125.0,0.6966));
  XSectionMap_wh.insert(std::pair<double,double>(125.5,0.6873));
  XSectionMap_wh.insert(std::pair<double,double>(126.0,0.6782));
  XSectionMap_wh.insert(std::pair<double,double>(126.5,0.6691));
  XSectionMap_wh.insert(std::pair<double,double>(127.0,0.6602));
  XSectionMap_wh.insert(std::pair<double,double>(127.5,0.6515));
  XSectionMap_wh.insert(std::pair<double,double>(128.0,0.6429));
  XSectionMap_wh.insert(std::pair<double,double>(128.5,0.6344));
  XSectionMap_wh.insert(std::pair<double,double>(129.0,0.6260));
  XSectionMap_wh.insert(std::pair<double,double>(129.5,0.6177));
  XSectionMap_wh.insert(std::pair<double,double>(130.0,0.6095));
  XSectionMap_wh.insert(std::pair<double,double>(130.5,0.6015));
  XSectionMap_wh.insert(std::pair<double,double>(131.0,0.5936));
  XSectionMap_wh.insert(std::pair<double,double>(131.5,0.5859));
  XSectionMap_wh.insert(std::pair<double,double>(132.0,0.5783));
  XSectionMap_wh.insert(std::pair<double,double>(132.5,0.5708));
  XSectionMap_wh.insert(std::pair<double,double>(133.0,0.5634));
  XSectionMap_wh.insert(std::pair<double,double>(133.5,0.5562));
  XSectionMap_wh.insert(std::pair<double,double>(134.0,0.5491));
  XSectionMap_wh.insert(std::pair<double,double>(134.5,0.5420));
  XSectionMap_wh.insert(std::pair<double,double>(135.0,0.5351));
  XSectionMap_wh.insert(std::pair<double,double>(135.5,0.5283));
  XSectionMap_wh.insert(std::pair<double,double>(136.0,0.5215));
  XSectionMap_wh.insert(std::pair<double,double>(136.5,0.5149));
  XSectionMap_wh.insert(std::pair<double,double>(137.0,0.5084));
  XSectionMap_wh.insert(std::pair<double,double>(137.5,0.5020));
  XSectionMap_wh.insert(std::pair<double,double>(138.0,0.4956));
  XSectionMap_wh.insert(std::pair<double,double>(138.5,0.4894));
  XSectionMap_wh.insert(std::pair<double,double>(139.0,0.4833));
  XSectionMap_wh.insert(std::pair<double,double>(139.5,0.4772));
  XSectionMap_wh.insert(std::pair<double,double>(140.0,0.4713));
  XSectionMap_wh.insert(std::pair<double,double>(141.0,0.4597));
  XSectionMap_wh.insert(std::pair<double,double>(142.0,0.4484));
  XSectionMap_wh.insert(std::pair<double,double>(143.0,0.4375));
  XSectionMap_wh.insert(std::pair<double,double>(144.0,0.4268));
  XSectionMap_wh.insert(std::pair<double,double>(145.0,0.4164));
  XSectionMap_wh.insert(std::pair<double,double>(146.0,0.4062));
  XSectionMap_wh.insert(std::pair<double,double>(147.0,0.3963));
  XSectionMap_wh.insert(std::pair<double,double>(148.0,0.3867));
  XSectionMap_wh.insert(std::pair<double,double>(149.0,0.3773));
  XSectionMap_wh.insert(std::pair<double,double>(150.0,0.3681));
    
  XSectionMap_zh.insert(std::pair<double,double>(110.0,0.5869));
  XSectionMap_zh.insert(std::pair<double,double>(110.5,0.5788));
  XSectionMap_zh.insert(std::pair<double,double>(111.0,0.5708));
  XSectionMap_zh.insert(std::pair<double,double>(111.5,0.5629));
  XSectionMap_zh.insert(std::pair<double,double>(112.0,0.5552));
  XSectionMap_zh.insert(std::pair<double,double>(112.5,0.5476));
  XSectionMap_zh.insert(std::pair<double,double>(113.0,0.5402));
  XSectionMap_zh.insert(std::pair<double,double>(113.5,0.5329));
  XSectionMap_zh.insert(std::pair<double,double>(114.0,0.5258));
  XSectionMap_zh.insert(std::pair<double,double>(114.5,0.5187));
  XSectionMap_zh.insert(std::pair<double,double>(115.0,0.5117));
  XSectionMap_zh.insert(std::pair<double,double>(115.5,0.5049));
  XSectionMap_zh.insert(std::pair<double,double>(116.0,0.4981));
  XSectionMap_zh.insert(std::pair<double,double>(116.5,0.4916));
  XSectionMap_zh.insert(std::pair<double,double>(117.0,0.4850));
  XSectionMap_zh.insert(std::pair<double,double>(117.5,0.4787));
  XSectionMap_zh.insert(std::pair<double,double>(118.0,0.4724));
  XSectionMap_zh.insert(std::pair<double,double>(118.5,0.4662));
  XSectionMap_zh.insert(std::pair<double,double>(119.0,0.4602));
  XSectionMap_zh.insert(std::pair<double,double>(119.5,0.4542));
  XSectionMap_zh.insert(std::pair<double,double>(120.0,0.4483));
  XSectionMap_zh.insert(std::pair<double,double>(120.5,0.4426));
  XSectionMap_zh.insert(std::pair<double,double>(121.0,0.4368));
  XSectionMap_zh.insert(std::pair<double,double>(121.5,0.4312));
  XSectionMap_zh.insert(std::pair<double,double>(122.0,0.4257));
  XSectionMap_zh.insert(std::pair<double,double>(122.5,0.4203));
  XSectionMap_zh.insert(std::pair<double,double>(123.0,0.4150));
  XSectionMap_zh.insert(std::pair<double,double>(123.5,0.4096));
  XSectionMap_zh.insert(std::pair<double,double>(124.0,0.4044));
  XSectionMap_zh.insert(std::pair<double,double>(124.5,0.3993));
  XSectionMap_zh.insert(std::pair<double,double>(125.0,0.3943));
  XSectionMap_zh.insert(std::pair<double,double>(125.5,0.3893));
  XSectionMap_zh.insert(std::pair<double,double>(126.0,0.3843));
  XSectionMap_zh.insert(std::pair<double,double>(126.5,0.3794));
  XSectionMap_zh.insert(std::pair<double,double>(127.0,0.3746));
  XSectionMap_zh.insert(std::pair<double,double>(127.5,0.3699));
  XSectionMap_zh.insert(std::pair<double,double>(128.0,0.3652));
  XSectionMap_zh.insert(std::pair<double,double>(128.5,0.3606));
  XSectionMap_zh.insert(std::pair<double,double>(129.0,0.3561));
  XSectionMap_zh.insert(std::pair<double,double>(129.5,0.3516));
  XSectionMap_zh.insert(std::pair<double,double>(130.0,0.3473));
  XSectionMap_zh.insert(std::pair<double,double>(130.5,0.3430));
  XSectionMap_zh.insert(std::pair<double,double>(131.0,0.3388));
  XSectionMap_zh.insert(std::pair<double,double>(131.5,0.3347));
  XSectionMap_zh.insert(std::pair<double,double>(132.0,0.3306));
  XSectionMap_zh.insert(std::pair<double,double>(132.5,0.3266));
  XSectionMap_zh.insert(std::pair<double,double>(133.0,0.3226));
  XSectionMap_zh.insert(std::pair<double,double>(133.5,0.3188));
  XSectionMap_zh.insert(std::pair<double,double>(134.0,0.3149));
  XSectionMap_zh.insert(std::pair<double,double>(134.5,0.3112));
  XSectionMap_zh.insert(std::pair<double,double>(135.0,0.3074));
  XSectionMap_zh.insert(std::pair<double,double>(135.5,0.3038));
  XSectionMap_zh.insert(std::pair<double,double>(136.0,0.3001));
  XSectionMap_zh.insert(std::pair<double,double>(136.5,0.2966));
  XSectionMap_zh.insert(std::pair<double,double>(137.0,0.2930));
  XSectionMap_zh.insert(std::pair<double,double>(137.5,0.2895));
  XSectionMap_zh.insert(std::pair<double,double>(138.0,0.2861));
  XSectionMap_zh.insert(std::pair<double,double>(138.5,0.2827));
  XSectionMap_zh.insert(std::pair<double,double>(139.0,0.2793));
  XSectionMap_zh.insert(std::pair<double,double>(139.5,0.2760));
  XSectionMap_zh.insert(std::pair<double,double>(140.0,0.2728));
  XSectionMap_zh.insert(std::pair<double,double>(141.0,0.2664));
  XSectionMap_zh.insert(std::pair<double,double>(142.0,0.2601));
  XSectionMap_zh.insert(std::pair<double,double>(143.0,0.2541));
  XSectionMap_zh.insert(std::pair<double,double>(144.0,0.2482));
  XSectionMap_zh.insert(std::pair<double,double>(145.0,0.2424));
  XSectionMap_zh.insert(std::pair<double,double>(146.0,0.2368));
  XSectionMap_zh.insert(std::pair<double,double>(147.0,0.2314));
  XSectionMap_zh.insert(std::pair<double,double>(148.0,0.2261));
  XSectionMap_zh.insert(std::pair<double,double>(149.0,0.2209));
  XSectionMap_zh.insert(std::pair<double,double>(150.0,0.2159));
    
  XSectionMap_tth.insert(std::pair<double,double>(110.0,0.1887));
  XSectionMap_tth.insert(std::pair<double,double>(110.5,0.1863));
  XSectionMap_tth.insert(std::pair<double,double>(111.0,0.1840));
  XSectionMap_tth.insert(std::pair<double,double>(111.5,0.1816));
  XSectionMap_tth.insert(std::pair<double,double>(112.0,0.1793));
  XSectionMap_tth.insert(std::pair<double,double>(112.5,0.1771));
  XSectionMap_tth.insert(std::pair<double,double>(113.0,0.1748));
  XSectionMap_tth.insert(std::pair<double,double>(113.5,0.1727));
  XSectionMap_tth.insert(std::pair<double,double>(114.0,0.1705));
  XSectionMap_tth.insert(std::pair<double,double>(114.5,0.1684));
  XSectionMap_tth.insert(std::pair<double,double>(115.0,0.1663));
  XSectionMap_tth.insert(std::pair<double,double>(115.5,0.1642));
  XSectionMap_tth.insert(std::pair<double,double>(116.0,0.1622));
  XSectionMap_tth.insert(std::pair<double,double>(116.5,0.1602));
  XSectionMap_tth.insert(std::pair<double,double>(117.0,0.1582));
  XSectionMap_tth.insert(std::pair<double,double>(117.5,0.1563));
  XSectionMap_tth.insert(std::pair<double,double>(118.0,0.1543));
  XSectionMap_tth.insert(std::pair<double,double>(118.5,0.1524));
  XSectionMap_tth.insert(std::pair<double,double>(119.0,0.1506));
  XSectionMap_tth.insert(std::pair<double,double>(119.5,0.1487));
  XSectionMap_tth.insert(std::pair<double,double>(120.0,0.1470));
  XSectionMap_tth.insert(std::pair<double,double>(120.5,0.1452));
  XSectionMap_tth.insert(std::pair<double,double>(121.0,0.1434));
  XSectionMap_tth.insert(std::pair<double,double>(121.5,0.1417));
  XSectionMap_tth.insert(std::pair<double,double>(122.0,0.1400));
  XSectionMap_tth.insert(std::pair<double,double>(122.5,0.1383));
  XSectionMap_tth.insert(std::pair<double,double>(123.0,0.1366));
  XSectionMap_tth.insert(std::pair<double,double>(123.5,0.1350));
  XSectionMap_tth.insert(std::pair<double,double>(124.0,0.1334));
  XSectionMap_tth.insert(std::pair<double,double>(124.5,0.1318));
  XSectionMap_tth.insert(std::pair<double,double>(125.0,0.1302));
  XSectionMap_tth.insert(std::pair<double,double>(125.5,0.1287));
  XSectionMap_tth.insert(std::pair<double,double>(126.0,0.1271));
  XSectionMap_tth.insert(std::pair<double,double>(126.5,0.1256));
  XSectionMap_tth.insert(std::pair<double,double>(127.0,0.1241));
  XSectionMap_tth.insert(std::pair<double,double>(127.5,0.1227));
  XSectionMap_tth.insert(std::pair<double,double>(128.0,0.1212));
  XSectionMap_tth.insert(std::pair<double,double>(128.5,0.1198));
  XSectionMap_tth.insert(std::pair<double,double>(129.0,0.1184));
  XSectionMap_tth.insert(std::pair<double,double>(129.5,0.1170));
  XSectionMap_tth.insert(std::pair<double,double>(130.0,0.1157));
  XSectionMap_tth.insert(std::pair<double,double>(130.5,0.1143));
  XSectionMap_tth.insert(std::pair<double,double>(131.0,0.1130));
  XSectionMap_tth.insert(std::pair<double,double>(131.5,0.1117));
  XSectionMap_tth.insert(std::pair<double,double>(132.0,0.1104));
  XSectionMap_tth.insert(std::pair<double,double>(132.5,0.1091));
  XSectionMap_tth.insert(std::pair<double,double>(133.0,0.1079));
  XSectionMap_tth.insert(std::pair<double,double>(133.5,0.1067));
  XSectionMap_tth.insert(std::pair<double,double>(134.0,0.1054));
  XSectionMap_tth.insert(std::pair<double,double>(134.5,0.1042));
  XSectionMap_tth.insert(std::pair<double,double>(135.0,0.1031));
  XSectionMap_tth.insert(std::pair<double,double>(135.5,0.1019));
  XSectionMap_tth.insert(std::pair<double,double>(136.0,0.1007));
  XSectionMap_tth.insert(std::pair<double,double>(136.5,0.09961));
  XSectionMap_tth.insert(std::pair<double,double>(137.0,0.09849));
  XSectionMap_tth.insert(std::pair<double,double>(137.5,0.09738));
  XSectionMap_tth.insert(std::pair<double,double>(138.0,0.09629));
  XSectionMap_tth.insert(std::pair<double,double>(138.5,0.09521));
  XSectionMap_tth.insert(std::pair<double,double>(139.0,0.09415));
  XSectionMap_tth.insert(std::pair<double,double>(139.5,0.09310));
  XSectionMap_tth.insert(std::pair<double,double>(140.0,0.09207));
  XSectionMap_tth.insert(std::pair<double,double>(141.0,0.09004));
  XSectionMap_tth.insert(std::pair<double,double>(142.0,0.08807));
  XSectionMap_tth.insert(std::pair<double,double>(143.0,0.08615));
  XSectionMap_tth.insert(std::pair<double,double>(144.0,0.08427));
  XSectionMap_tth.insert(std::pair<double,double>(145.0,0.08246));
  XSectionMap_tth.insert(std::pair<double,double>(146.0,0.08068));
  XSectionMap_tth.insert(std::pair<double,double>(147.0,0.07895));
  XSectionMap_tth.insert(std::pair<double,double>(148.0,0.07727));
  XSectionMap_tth.insert(std::pair<double,double>(149.0,0.07563));
  XSectionMap_tth.insert(std::pair<double,double>(150.0,0.07403));
    
    
  for(int i=0; i<numsmpoints; ++i) {
    double mass = 110. + i*0.5;
    if (mass < 140. || mass == (double)(int)mass) {

/*       std::cout<<" assignig CS for mass "<<mass<<" normaly  as "<<XSectionMap_ggh[mass]<<std::endl; */
/*       std::cout<<"                                             "<<XSectionMap_vbf[mass]<<std::endl; */
/*       std::cout<<"                                             "<<XSectionMap_wh[mass]<<std::endl; */
/*       std::cout<<"                                             "<<XSectionMap_zh[mass]<<std::endl; */
/*       std::cout<<"                                             "<<XSectionMap_tth[mass]<<std::endl; */

      gghxsec[i] = XSectionMap_ggh[mass];
      vbfxsec[i] = XSectionMap_vbf[mass];
      whxsec [i] = XSectionMap_wh[mass];
      zhxsec [i] = XSectionMap_zh[mass];
      tthxsec[i] = XSectionMap_tth[mass];
    } else {

/*       std::cout<<"  injterpolated CS for mass = "<<mass<<" is "<<(XSectionMap_ggh[mass-0.5]+XSectionMap_ggh[mass+0.5])/2.<<std::endl; */
/*       std::cout<<"                                            "<<(XSectionMap_vbf[mass-0.5]+XSectionMap_vbf[mass+0.5])/2.<<std::endl; */
/*       std::cout<<"                                            "<<(XSectionMap_wh [mass-0.5]+XSectionMap_wh [mass+0.5])/2.<<std::endl; */
/*       std::cout<<"                                            "<<(XSectionMap_zh [mass-0.5]+XSectionMap_zh [mass+0.5])/2.<<std::endl; */
/*       std::cout<<"                                            "<<(XSectionMap_tth[mass-0.5]+XSectionMap_tth[mass+0.5])/2.<<std::endl; */

      gghxsec[i] = (XSectionMap_ggh[mass-0.5]+XSectionMap_ggh[mass+0.5])/2.;
      vbfxsec[i] = (XSectionMap_vbf[mass-0.5]+XSectionMap_vbf[mass+0.5])/2.;
      whxsec [i] = (XSectionMap_wh [mass-0.5]+XSectionMap_wh [mass+0.5])/2.;
      zhxsec [i] = (XSectionMap_zh [mass-0.5]+XSectionMap_zh [mass+0.5])/2.;
      tthxsec[i] = (XSectionMap_tth[mass-0.5]+XSectionMap_tth[mass+0.5])/2.;
    }
  }
    
    
  for(int i=0; i<numsmpoints; ++i) {
    wzhxsec[i] = (whxsec[i] + zhxsec[i]);
    //tthxsec[i] = tthxsec[i];

    // interference with 2photon BG, -2.5% on CS in all masses on GF process
    gghxsec[i] = 0.975 * gghxsec[i];
  }

  processCrossSectionMap.insert(std::pair<TString,double*>("gghxsec",gghxsec));
  processCrossSectionMap.insert(std::pair<TString,double*>("sm4gghxsec",sm4gghxsec));
  processCrossSectionMap.insert(std::pair<TString,double*>("vbfxsec",vbfxsec));
  
  processCrossSectionMap.insert(std::pair<TString,double*>("wzhxsec",wzhxsec));
  
  processCrossSectionMap.insert(std::pair<TString,double*>("tthxsec",tthxsec));
  
  processCrossSectionMap.insert(std::pair<TString,double*>("smbr",smbr));
  processCrossSectionMap.insert(std::pair<TString,double*>("ffbr",ffbr));
  processCrossSectionMap.insert(std::pair<TString,double*>("sm4br",sm4br));

  for(int iMass=0; iMass<numsmpoints; ++iMass) nullarray[iMass] = 0.;
  processCrossSectionMap.insert(std::pair<TString,double*>("null",nullarray));
  
  return;
}

#endif
