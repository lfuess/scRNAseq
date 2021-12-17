##here we will conduct the stats shown in supplemental table 3 in our MS
##these stats test for enrichment of immune cell markers general, and type-specific enrichment
##we look at three different MS: Fuess et al. 2021a & b, and Lohman 2017
##and we use the following formula to test enrichment + a prop test for ones where sig enrichment is detected
  #dattable <- matrix(c(markers + sig tau, total sig tau, total markers, total genes), nrow = 2)


##first we'll look for over-representation of immune cell markers generally, and specific immune cell types
##in our gene expression-microbiome association MS (Fuess 2021; mbio)##


  ##first we look and see if immune cell type marker genes are over-represented generally##
  ##beta##
  dattable = matrix(c(20,289,639,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##caulo##
  dattable = matrix(c(30,289,591,17231),nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##chitin##
  dattable = matrix(c(17,289,389,17231),nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##chlam##
  dattable = matrix(c(20,289,376,17231),nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##clostrid##
  dattable = matrix(c(31,289,926,17231),nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##diversity##
  dattable = matrix(c(64,289,1928,17231),nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##geo##
  dattable = matrix(c(18,289,661,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##gp##
  dattable = matrix(c(11,289,505,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##halo##
  dattable = matrix(c(17,289,393,17231), nrow=2)
  chisq.test(dattable)
  ##sig
  
  #inc##
  dattable = matrix(c(23,289,480,17231), nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##methyl##
  dattable = matrix(c(22,289,523,17231), nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##nocard##
  dattable = matrix(c(33,289,892,17231), nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##orb##
  dattable = matrix(c(11,289,306,17231), nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##pepto##
  dattable = matrix(c(33,289,908,17231), nrow=2)
  chisq.test(dattable)
  ##sig
  
  ##sparto##
  dattable = matrix(c(14,289,488,17231), nrow=2)
  chisq.test(dattable)
  ##not sig
  
  ##rubro##
  dattable = matrix(c(12,289,351,17231), nrow=2)
  chisq.test(dattable)
  ##sig


  ##then for groups with significant over-representation of immune cell markers generally, we look for over-representation of specific cell types##
  ##Start with beta##
  ##APC##
  dattable = matrix(c(1,59,639,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(3,38,639,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(2,37,639,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##hc##
  dattable = matrix(c(1,29,639,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(12,95,639,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##RBC##
  dattable = matrix(c(1,18,639,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##Diversity##
  ##APC##
  dattable = matrix(c(11,59,1928,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(16,38,1928,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(14,37,1928,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(5,29,1928,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(14,95,1928,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(1,8,1928,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##platelet##
  dattable = matrix(c(2,16,1928,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(1,18,1928,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##Caulo##
  ##APC##
  dattable = matrix(c(4,59,591,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(11,38,591,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(6,37,591,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(2,29,591,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(6,95,591,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(1,8,591,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##Chitin##
  ##APC##
  dattable = matrix(c(2,59,389,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(4,38,389,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(8,37,389,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(2,29,389,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(1,95,389,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##Chlam##
  ##APC##
  dattable = matrix(c(1,59,376,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(4,38,376,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(5,37,376,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##neut##
  dattable = matrix(c(8,95,376,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##NKC##
  dattable = matrix(c(1,8,376,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(1,18,376,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##Clostrid##
  ##APC##
  dattable = matrix(c(2,59,926,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(5,38,926,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(6,37,926,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(1,29,926,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(16,95,926,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##NKC##
  dattable = matrix(c(1,8,926,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##Geo##
  ##APC##
  dattable = matrix(c(2,59,661,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(4,37,661,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(11,95,661,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##RBC##
  dattable = matrix(c(1,18,661,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##GP10##
  ##APC##
  dattable = matrix(c(3,59,505,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(2,37,505,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##hc##
  dattable = matrix(c(2,29,505,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(4,95,505,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##Halo##
  ##APC##
  dattable = matrix(c(1,59,393,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(3,38,393,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(6,37,393,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(1,29,393,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(6,95,393,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  
  ##Inc##
  ##APC##
  dattable = matrix(c(1,59,480,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(5,38,480,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(3,37,480,17231),nrow=2)
  chisq.test(dattable)
  ##notsig##
  
  ##hc##
  dattable = matrix(c(1,29,480,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(11,95,480,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##NKC##
  dattable = matrix(c(1,8,480,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(1,18,480,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##Methyl##
  ##APC##
  dattable = matrix(c(2,59,523,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(5,38,523,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(8,37,523,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(1,29,523,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(6,95,523,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  
  ##Nocard##
  ##APC##
  dattable = matrix(c(2,59,892,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(7,38,892,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(7,37,892,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(1,29,892,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(14,95,892,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##NKC##
  dattable = matrix(c(1,8,892,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(1,18,892,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  
  ##Orb##
  ##APC##
  dattable = matrix(c(1,59,306,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(2,38,306,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##hc##
  dattable = matrix(c(1,29,306,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(7,95,306,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  
  ##Pepto##
  ##APC##
  dattable = matrix(c(3,59,908,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(7,38,908,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(6,37,908,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(2,29,908,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(14,95,908,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##NKC##
  dattable = matrix(c(1,8,908,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##Rubro##
  ##APC##
  dattable = matrix(c(1,59,351,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(1,38,351,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(4,37,351,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##neut##
  dattable = matrix(c(6,95,351,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  
  ##Sparto##
  
  ##B-cell##
  dattable = matrix(c(2,38,488,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(5,37,488,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(2,29,488,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(5,95,488,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##Prop tests for the above sig enrichment groups##
    ##mbio-all diversity
    prop.test(23,64)
    ##mbio-Bcell diversity
    prop.test(4,16)
    ##mbio-fibroblasts diversity
    prop.test(3,14)
    ##mbio-all beta
    prop.test(2,20)
    ##mbio-neutrophil beta
    prop.test(12,12)
    ##mbio-all caulo
    prop.test(4,30)
    ##mbio-bcell caulo
    prop.test(1,11)
    ##mbio-fibroblast caulo
    prop.test(2,6)
    ##mbio-all chitin
    prop.test(4,17)
    ##mbio-Bcell chitin
    prop.test(4,4)
    ##mbio-fibro chitin
    prop.test(3,8)
    ##mbio-all chlam
    prop.test(5,20)
    ##mbio-neut chlam
    prop.test(8,8)
    ##mbio-bcells chlam
    prop.test(3,4)
    ##mbio-fib chlam
    prop.test(4,5)
    ##mbio-all clostrid
    prop.test(25,31)
    ##mbio-neut clostrid
    prop.test(16,16)
    ##mbio-fib clostrid
    prop.test(5,6)
    ##mbio-neut geo
    prop.test(11,11)
    ##mbio-all halo
    prop.test(17,18)  
    ##mbio-neut halo
    prop.test(6,6)  
    ##mbio-fib halo
    prop.test(5,6) 
    ##mbio-all inc
    prop.test(7,23) 
    ##mbio-neut inc
    prop.test(11,11) 
    ##mbio-bcell inc
    prop.test(4,5) 
    ##mbio-all nocard
    prop.test(23,33) 
    ##mbio-neut nocard
    prop.test(14, 14) 
    ##mbio-bcells nocard
    prop.test(5, 7) 
    ##mbio-all methyl
    prop.test(4, 22) 
    ##mbio-bcells methyl
    prop.test(5, 5) 
    ##mbio-fib methyl
    prop.test(3, 8) 
    ##mbio-all orb
    prop.test(11, 11) 
    ##mbio-neut orb
    prop.test(7, 7) 
    ##mbio-all pepto
    prop.test(23, 33) 
    ##mbio-neut pepto
    prop.test(14, 14) 
    ##mbio-bcell pepto
    prop.test(6, 7) 
    ##mbio-fib pepto
    prop.test(5, 6) 
    ##mbio-all rubro
    prop.test(12, 12) 
    ##mbio-neut rubro
    prop.test(6, 6) 
    ##mbio-fib spart
    prop.test(2, 5) 



##Next we'll look for significant over-representation of markers (generally) and
##cell-type specific marerks in the Fuess et al. 2021 TagSeq MS (molec ecol)
  ##infection-general##
  dattable = matrix(c(105,289,2369,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##APC##
  dattable = matrix(c(37,59,2369,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##B-cell##
  dattable = matrix(c(14,38,2369,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(12,37,2369,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(12,29,2369,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##neut##
  dattable = matrix(c(20,95,2369,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(4,8,2369,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##platelet##
  dattable = matrix(c(5,16,2369,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(1,18,2369,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  
  ##fibrosis-nothing is significant!!!##
  dattable = matrix(c(119,289,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##APC##
  dattable = matrix(c(23,59,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(14,38,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(18,37,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##hc##
  dattable = matrix(c(4,29,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##neut##
  dattable = matrix(c(43,95,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(5,8,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##platelet##
  dattable = matrix(c(6,16,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(6,18,5826,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  
  ##RBCvGBC##
  dattable = matrix(c(44,289,1202,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##APC##
  dattable = matrix(c(5,59,1202,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(15,38,1202,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##fibroblast##
  dattable = matrix(c(2,37,1202,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##hc##
  dattable = matrix(c(8,29,1202,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##neut##
  dattable = matrix(c(12,95,1202,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##platelet##
  dattable = matrix(c(1,16,1202,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(1,18,1202,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  
  ##GvF##
  dattable = matrix(c(135,289,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##APC##
  dattable = matrix(c(21,59,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(15,38,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(23,37,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##hc##
  dattable = matrix(c(8,29,7745,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##neut##
  dattable = matrix(c(47,95,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(2,8,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##platelet##
  dattable = matrix(c(9,16,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(10,18,7745,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  
  ##FvR##
  dattable = matrix(c(223,289,10601,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##APC##
  dattable = matrix(c(28,59,10601,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(34,38,10601,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(26,37,10601,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##hc##
  dattable = matrix(c(24,29,10601,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##neut##
  dattable = matrix(c(82,95,10601,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##NKC##
  dattable = matrix(c(3,8,10601,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##platelet##
  dattable = matrix(c(12,16,10601,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##RBC##
  dattable = matrix(c(14,18,10601,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  
  ##GvRxI##
  dattable = matrix(c(1,289,4,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(1,8,4,17231),nrow=2)
  chisq.test(dattable)
  ##sig...technically##
  
  ##FvGxI##
  dattable = matrix(c(14,289,558,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##APC##
  dattable = matrix(c(3,59,558,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(3,38,558,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(2,37,558,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##hc##
  dattable = matrix(c(1,29,558,17231),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  ##neut##
  dattable = matrix(c(2,95,558,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(2,8,558,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##platelet##
  dattable = matrix(c(1,16,558,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##PopxFib##
  dattable = matrix(c(5,289,84,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##APC##
  dattable = matrix(c(2,59,84,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##B-cell##
  dattable = matrix(c(1,38,84,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(1,37,84,17231),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##NKC##
  dattable = matrix(c(1,8,84,17231),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##next we need to check for significant prop differences##
    ##fuess- all markers infection
    prop.test(44,105)
    ##fuess- HCs infection
    prop.test(3,12)
    ##fuess- APCs infection
    prop.test(37,37)
    ##fuess- B-cells infection
    prop.test(13,14)
    ##fuess- fibroblasts infection
    prop.test(7,12)
    ##fuess-all markers RvG
    prop.test(5,44)
    ##fuess-HCs RvG
    prop.test(8,8)
    ##fuess-Bcells RvG
    prop.test(15,15)
    ##fuess-all markers FvR
    prop.test(178,223)
    ##fuess-Neut FvR
    prop.test(75,82)
    ##fuess-all markers Pop*fib
    prop.test(4,5)
    ##fuess-APC Pop*fib
    prop.test(2,2)
    ##fuess-NKC Pop*fib
    prop.test(1,1)

##finally we look at enrichment in Lohman 2017 data set##
  ##infection##
  dattable = matrix(c(10,289,64,9078),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##APC##
  dattable = matrix(c(2,59,64,9078),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##RBC##
  dattable = matrix(c(5,18,64,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##platelet##
  dattable = matrix(c(1,16,64,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(1,37,64,9078),nrow=2)
  chisq.test(dattable)
  ##barely not sig##
  
  
  ##popuation##
  dattable = matrix(c(43,289,643,9078),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##hc##
  dattable = matrix(c(10,29,643,9078),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##neut##
  dattable = matrix(c(14,95,643,9078),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##APC##
  dattable = matrix(c(8,59,643,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##B-cell##
  dattable = matrix(c(6,38,643,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##platelet##
  dattable = matrix(c(1,16,643,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##fibroblast##
  dattable = matrix(c(2,37,643,9078),nrow=2)
  chisq.test(dattable)
  ## not sig##
  
  ##NKC##
  dattable = matrix(c(2,8,643,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##
  
  ##interaction##
  dattable = matrix(c(1,289,16,9078),nrow=2)
  chisq.test(dattable)
  ##sig##
  
  ##APC##
  dattable = matrix(c(1,59,16,9078),nrow=2)
  chisq.test(dattable)
  ##not sig##

  ##finally check for directional enrichment##
    ##lohman-all markers infection##
    prop.test(1, 10)
    ##lohman- RBC markers, infection
    prop.test(0,6)
    ##lohman-all markers population##
    prop.test(20, 43)
    ##lohman- neut markers, population
    prop.test(5,14)
    ##lohman- HC markers, population
    prop.test(1,10)

