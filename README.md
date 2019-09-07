
Cascade
-------
Software utility to manipulate touchstone s-parameter files.

Features
--------
The utility takes a touchstone file as standard input
and outputs a new touchstone file (or other formats)
to standard output.  Cascade provides the following
input transformations as command line options:

```
-cascade             : cascade together the two networks on top of stack
-deembed             : de-embed the output of the top two networks on the stack
-ideembed            : de-embed the input of the top two networks on the stack
-swap                : swap the top two networks on the stack
-cbg                 : transform network on top of stack into a common-base arrangement
-ccd                 : transform network on top of stack into a common-collector arrangement
-lift <a+bj>         : lift network on top of stack from ground and insert an impedance, j required
-lift <henries>      : lift network on top of stack from ground and insert an inductor
-p                   : print network on top of stack
-series <ohms>       : push a series resistor onto the stack
-shunt <ohms>        : push a shunt resistor onto the stack
-f <filename>        : push a touchstone network onto the stack
```

By default the utility writes out the network on the top of
the stack in touchstone format with GUM and Rollet stability information 
as comments.  It can also output the network in following alternative formats:

```
-n  : print the result of write_touchstone() from scikit-rf 
-a  : display the network as ABCD matrices
-s  : summarize the network in terms of impedances, stability and gain (dB) values
```

Only 50 ohm networks are supported.

Installation
------------
```
$ sudo pip install -r requirements.txt
$ sudo cp cascade.py /usr/local/bin/cascade
$ sudo chmod 755 /usr/local/bin/cascade
$ sh test.sh       # to run unit tests
```

Examples
--------

Add GUM (dB), K and mu as comments to touchstone file.


```
$ < 2n5179_5ma.s2p cascade
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.471    0.00      6.78  180.00         0    0.00     0.844    0.00 !  23.1       inf     1.185
100       0.471  -90.00      6.78  122.00     0.023   64.00     0.844  -51.00 !  23.1    0.4623    0.8889
200       0.314 -145.00       4.2  100.00     0.034   58.00      0.78  -93.00 !  17.0     1.109     1.021
300        0.23  156.00      2.76   91.00     0.043   65.00     0.768 -134.00 !  12.9     1.819     1.145
400       0.171  108.00      2.17   86.00     0.056   63.00     0.756 -177.00 !  10.5     1.874     1.157
500       0.168   54.00      1.86   79.00     0.062   62.00     0.741  140.00 !   9.0     1.883     1.147
600       0.149   -9.00      1.53   71.00     0.069   66.00      0.74   98.00 !   7.2     2.074     1.164
700       0.137  -72.00      1.31   67.00     0.073   71.00     0.739   54.00 !   5.9     2.469     1.213
800       0.119 -129.00      1.18   64.00     0.092   74.00     0.744    8.00 !   5.0     2.098     1.174
900       0.153 -174.00      1.13   58.00     0.101   68.00     0.742  -38.00 !   4.6     1.875     1.142
1000      0.171  122.00     0.979   49.00     0.106   71.00     0.749  -82.00 !   3.5     2.083     1.164
```


Display the result as ABCD matrices.


```
$ < 2n5179_5ma.s2p cascade -a
MHZ            A                 B                 C                 D
0       0.01692 -180.00        10 -180.00 0.0001217 -180.00   0.07194 -180.00
100     0.05534  -88.08     7.139 -166.65  0.001397  -51.43    0.1243 -120.34
200       0.107  -70.35     6.367 -148.80  0.004076  -59.44      0.17 -131.82
300      0.2218  -60.14     6.235 -130.84  0.007505  -77.69    0.1408 -149.65
400      0.3684  -72.05     4.008  -95.97  0.009134  -94.89   0.03424  -93.31
500      0.4608  -86.97     9.595  -33.46  0.008633 -105.44    0.1543  -28.46
600      0.4624 -105.15     22.05  -37.35  0.008063 -103.97    0.3211  -24.06
700      0.2943 -122.15     32.17  -54.24  0.006801 -106.18    0.5633  -33.08
800     0.06799  -81.22     36.07  -68.63   0.00344  -87.15    0.7668  -53.02
900       0.246    0.28        33  -77.79  0.006156  -18.70    0.8021  -71.03
1000     0.5662    4.77     33.48  -76.44   0.01307  -21.60    0.6974  -89.08
```


Summarize the network.


```
$ < 2n5179_5ma.s2p cascade -s
MHZ            ZIN             ZOUT        GUM    GUI    S21    GUO   GMSG   GMAG     GU         K        MU
0               139+0j           591+0j  23.13   1.09  16.62   5.41    inf  23.13   0.00       inf     1.185
100       31.84-38.55j     22.13-100.9j  23.13   1.09  16.62   5.41  24.70      -      -    0.4623    0.8889
200       27.94-11.17j     11.59-46.09j  16.99   0.45  12.46   4.07  20.92  18.91   0.71     1.109     1.021
300        32.15+6.35j     7.719-20.79j  12.92   0.24   8.82   3.87  18.07  12.84  -0.33     1.819     1.145
400       42.77+14.33j     6.952-1.284j  10.54   0.13   6.73   3.68  15.88  10.49  -0.26     1.874     1.157
500       58.49+16.36j     8.399+17.74j   8.97   0.12   5.39   3.46  14.77   9.36   0.17     1.883     1.147
600       67.17-3.202j      12.9+41.79j   7.24   0.10   3.69   3.44  13.46   7.56   0.15     2.074     1.164
700       52.52-13.95j      33.5+88.26j   5.86   0.08   2.35   3.43  12.54   5.79  -0.17     2.469     1.213
800       42.35-7.945j       279+129.4j   5.00   0.06   1.44   3.50  11.08   5.12  -0.04     2.098     1.174
900       36.78-1.205j     58.96-119.9j   4.64   0.10   1.06   3.47  10.49   5.10   0.24     1.875     1.142
1000       40.1+11.98j     16.23-54.84j   3.52   0.13  -0.18   3.58   9.65   3.73   0.04     2.083     1.164
```


Shunt a 330 ohm resistor across the output of the two-port network.  This makes the transistor unconditionally stable.


```
$ < 2n5179_5ma.s2p cascade -shunt 330 -cascade
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.471    0.00     5.949  180.00         0    0.00     0.618    0.00 !  18.7       inf     1.618
100      0.4695  -88.72     6.069  124.55   0.02059   66.55    0.6577  -53.05 !  19.2     1.557     1.129
200      0.3082 -143.49      3.91  103.15   0.03165   61.15    0.6784  -95.81 !  15.0     1.884     1.182
300      0.2213  155.91     2.664   93.31    0.0415   67.31    0.7377 -135.77 !  12.1     2.142     1.199
400      0.1643  105.92      2.13   86.17   0.05498   63.17    0.7603 -177.12 !  10.4     1.906     1.155
500      0.1675   51.12       1.8   77.00      0.06   60.00    0.7204  141.75 !   8.4     2.151     1.187
600       0.155  -10.66     1.431   68.02   0.06452   63.02    0.6532  101.19 !   5.6         3     1.318
700      0.1428  -70.78     1.181   64.66   0.06579   68.66    0.5781   57.14 !   3.3     4.307     1.543
800      0.1197 -125.52     1.043   63.60    0.0813   73.60    0.5415    8.59 !   1.9     4.169     1.601
900      0.1491 -171.48     1.008   59.77   0.09013   69.77    0.5603  -40.49 !   1.8     3.657     1.501
1000     0.1638  121.96    0.9022   51.97   0.09769   73.97     0.632  -85.31 !   1.4     3.358     1.383
```


Cascade a series 20 ohm resistor with the output of the two-port network.


```
$ < 2n5179_5ma.s2p cascade -series 20 -cascade
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.471    0.00     6.575  180.00         0    0.00    0.8487    0.00 !  23.0       inf     1.178
100      0.4714  -93.44     6.155  115.16   0.02088   57.16    0.7407  -46.86 !  20.3     1.194     1.041
200      0.3248 -148.73     3.448   92.65   0.02791   50.65    0.5297  -82.32 !  12.7     3.225      1.49
300       0.248  155.65     2.105   86.17   0.03279   60.17    0.3788 -122.00 !   7.4     5.894     2.244
400       0.185  111.63     1.606   85.66   0.04145   62.66    0.2999 -175.86 !   4.7     6.681     2.792
500      0.1686   59.96     1.412   83.15   0.04708   66.15    0.3478  127.85 !   3.7      6.41     2.351
600      0.1339   -5.65     1.245   77.85   0.05613   72.85    0.4864   85.50 !   3.2      5.38     1.783
700      0.1236  -76.67      1.17   73.13   0.06521   77.13    0.6494   47.27 !   3.8      3.85     1.394
800      0.1201 -138.88     1.121   65.13   0.08738   75.13    0.7539    7.12 !   4.7     2.218     1.171
900       0.163  179.31      1.04   53.18   0.09292   63.18    0.7004  -33.51 !   3.4     2.516     1.229
1000     0.1883  121.15    0.8238   41.83   0.08919   63.83    0.5533  -71.64 !   0.1     4.533     1.568
```


Lift terminal 3, the port connected to ground, and add a 10nH inductor.  Then cascade the result
with a shunt 100 ohm resistor to stabilize the result.


```
$ < 2n5179_5ma.s2p cascade -lift 10e-9 -shunt 100 -cascade
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.471    0.00     4.641  180.00 1.187e-18 -180.00    0.2621   -0.00 !  14.7 6.576e+16     3.815
100      0.2894  -31.98      3.67  112.81   0.03295  129.17    0.2602  -55.49 !  12.0     3.527     2.349
200      0.2813    0.15     2.299   91.42    0.1193  112.36    0.2716  -97.16 !   7.9     1.653     1.612
300      0.3553    4.07     1.669   75.50    0.2257   82.69    0.3018 -146.61 !   5.4     1.188       1.2
400      0.3673   -2.21     1.266   61.02    0.2619   55.24    0.2773  162.64 !   3.0     1.297     1.292
500      0.3933   -3.40    0.9646   51.78    0.2325   39.45    0.2338  118.67 !   0.7      1.81     1.675
600      0.3965  -14.35    0.7704   47.45     0.207   37.88    0.2079   80.31 !  -1.3      2.53     2.132
700       0.319  -21.21    0.6623   46.78    0.1807   37.44    0.2023   41.51 !  -2.9     3.631     2.797
800      0.2797   -3.83    0.5861   50.23    0.1252   54.88    0.2335    0.19 !  -4.0     5.997     3.258
900       0.375    0.32    0.6423   57.44    0.2209  100.88    0.3044  -46.80 !  -2.8     2.922     2.412
1000     0.4493  -17.00    0.8127   50.80    0.4842   83.82    0.3872 -103.62 !  -0.1     1.103     1.127
```


Insert a one ohm resistor at the emitter to provide shunt feedback.


```
$ < 2n5179_5ma.s2p cascade -lift 1+0j
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.507   -0.00     6.308  180.00 0.0007679    0.00    0.8541   -0.00 !  23.0     21.18     1.166
100      0.5026  -82.33     6.453  124.47   0.03084   68.45    0.8719  -48.35 !  23.7    0.1972    0.8006
200      0.3182 -135.30     4.084  101.58   0.04909   53.96    0.7981  -89.87 !  17.1    0.6767    0.9192
300       0.207  163.41     2.689   91.88   0.05657   50.07    0.7621 -130.87 !  12.6     1.447     1.098
400      0.1457  109.26     2.112   86.78   0.06327   47.44    0.7312 -174.12 !   9.9     1.872     1.178
500      0.1549   49.71     1.817   79.92   0.06313   48.98    0.7097  142.15 !   8.3     2.119     1.197
600      0.1538  -13.14     1.505   71.80    0.0672   56.96    0.7137   99.11 !   6.7     2.338     1.206
700      0.1477  -71.33     1.295   67.61   0.06973   65.17     0.724   54.26 !   5.6     2.716     1.243
800      0.1235 -123.93     1.171   64.46   0.09085   72.14    0.7409    7.85 !   4.9     2.168     1.183
900      0.1491 -169.65     1.126   58.22    0.1062   66.54    0.7445  -38.10 !   4.6     1.784     1.133
1000     0.1609  123.17    0.9781   48.86    0.1141   66.76    0.7477  -81.82 !   3.5     1.947     1.153
```


Transform the common emitter s-parameter input into a common-base.


```
$ < 2n5179_5ma.s2p cascade -cbg | tee cb.s2p
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0        0.6692  180.00      1.55   -0.00    0.0258    0.00    0.9649   -0.00 !  18.0     1.141     1.011
100      0.7455  173.47     1.673  -11.49   0.08774   75.40     1.046  -10.87 !   nan   0.04985    0.7311
200      0.7834  166.09     1.623  -22.23    0.1673   79.26     1.115  -27.69 !   nan   -0.1685    0.4953
300       0.798  154.26      1.41  -32.60    0.3239   72.65     1.017  -53.92 !   nan   -0.1879    0.3329
400      0.7442  144.59     1.256  -34.55     0.464   48.66    0.6954  -75.56 !   8.4    0.1916    0.3582
500      0.7206  135.76     1.264  -36.70     0.528   27.09    0.3917  -83.59 !   5.9    0.5052    0.3812
600      0.6681  125.60     1.315  -44.49    0.5151    3.26    0.1657  -65.59 !   5.1    0.7577    0.4936
700      0.5996  116.64     1.338  -57.47    0.3734  -19.62    0.2935    7.81 !   4.9    0.9993    0.9977
800      0.5371  109.65     1.212  -74.72    0.1278  -28.04    0.7227   -6.39 !   6.4     1.523     1.194
900      0.4694  105.01    0.9932  -93.71     0.209   81.20     1.083  -40.18 !   nan   -0.3312    0.7231
1000      0.414   89.52    0.5058 -123.81    0.5353   63.17     1.108  -80.94 !   nan   -0.3771    0.6711
```


Create a cascode amplifier.  S12 is significantly reduced compared to a CE amp.


```
$ < 2n5179_5ma.s2p cascade -f cb.s2p -cascade
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.471    0.00     6.716  180.00         0    0.00    0.9865    0.00 !  33.3       inf     1.014
100      0.4483  -80.32     7.881  132.15  0.001402  161.04     1.108   -7.70 !   nan    -8.698    0.8867
200      0.2328 -128.95     6.756  113.18  0.005638  172.67     1.305  -23.48 !   nan    -8.781    0.7421
300      0.0308  150.86     8.191   84.93   0.02932  164.18     1.743  -59.59 !   nan    -4.073    0.5058
400     0.07062  168.43     4.502   21.58   0.04292   81.79    0.7419 -136.29 !  16.6     1.237     1.063
500     0.09522   43.22     2.166   12.99   0.03016   59.78   0.07264  127.14 !   6.8     7.577     7.245
600      0.1466  -28.60     1.437   12.42   0.02538   55.17    0.3441   15.47 !   3.8     11.84     2.616
700      0.1725  -78.73     1.218   12.40   0.01894   54.25    0.5338   -5.24 !   3.3     15.13     1.805
800      0.1605 -119.73     1.156    5.91  0.009504   62.59    0.7571  -13.10 !   5.1     19.04     1.304
900      0.1849 -157.40     1.218  -15.35   0.02291  169.56     1.248  -38.84 !   nan    -9.861      0.78
1000      0.181  142.01    0.7138  -71.46   0.08179  137.52     1.286  -92.10 !   nan    -5.604    0.7372
```


Stabilize the cascode amp with a 100 ohm resistor across the output.


```
$ < 2n5179_5ma.s2p cascade -f cb.s2p -cascade -shunt 100 -cascade
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
0         0.471    0.00     4.487  180.00         0    0.00    0.3273    0.00 !  14.6       inf     3.055
100      0.4465  -80.38     5.168  133.55 0.0009194  162.44    0.3832   -9.59 !  15.9     71.73     2.554
200      0.2298 -130.26     4.346  117.98  0.003626  177.47    0.4864  -26.26 !  14.2     22.95     1.988
300     0.05872  112.41     5.397   99.27   0.01931  178.52    0.7909  -55.60 !  18.9     1.877     1.122
400     0.06053 -154.28     4.008   28.13   0.03821   88.34    0.6741 -142.93 !  14.7     1.828       1.2
500     0.08392   38.87     1.748   12.32   0.02434   59.11     0.231  170.60 !   5.1     11.08     3.688
600      0.1474  -31.25     1.078   11.43   0.01904   54.18   0.05164   90.05 !   0.8     23.77     13.21
700       0.176  -79.49    0.8808   12.91    0.0137   54.76    0.1107  -13.32 !  -0.9     39.69     8.258
800      0.1624 -119.61    0.8056    7.62  0.006623   64.30    0.2294  -21.30 !  -1.5     86.48     4.271
900      0.1815 -156.44    0.8088   -7.88   0.01521  177.03    0.5032  -43.32 !  -0.4     29.27      1.93
1000     0.1758  145.27     0.558  -56.91   0.06394  152.07    0.7856  -91.95 !  -0.8     5.114     1.209
```


