#RNA-seq Cgigas


For this project I will try to find differentially expressed genes in oysters (3) exposed to heat.

Overall steps - will use DESeq2-DAVID-REvigo  and SQLSHARE.


#Count Data

Count data was from CLC and is located in the data directory.

```
Stevens-MacBook-Air-3:Desktop sr320$ cd ASI-rna-seq/
Stevens-MacBook-Air-3:ASI-rna-seq sr320$ ls
data		notebook	output		scripts
Stevens-MacBook-Air-3:ASI-rna-seq sr320$ cd data
Stevens-MacBook-Air-3:data sr320$ curl -O -k https://raw.githubusercontent.com/sr320/austral/master/modules/data/Cgigas-HS-count.txt
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  863k  100  863k    0     0   304k      0  0:00:02  0:00:02 --:--:--  357k
Stevens-MacBook-Air-3:data sr320$ 
```

```
$ head Cgigas-HS-count.txt 
Feature ID	2M-Total	4M-Total	6m-Total	2M-HS-Total	4M-HS-Total	6M-HS-Total
CGI_10000001	25	26	35	56	10	35
CGI_10000002	56	16	36	34	19	29
CGI_10000003	0	19	25	6	79	8
CGI_10000004	0	2	1	5	0	2
CGI_10000005	0	0	0	0	0	0
CGI_10000009	93	68	58	1384	1287	816
CGI_10000010	185	195	254	762	774	1326
CGI_10000011	31	2	10	131	19	17
CGI_10000012	0	0	0	0	0	0
```




#Run DESeq in R


I will be using this page as a reference for analysis.

<http://nbviewer.ipython.org/github/sr320/eimd-sswd/blob/master/eimd_analysis.ipynb>

If running for first time - 

```
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

```