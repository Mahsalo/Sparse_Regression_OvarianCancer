# Ovarian-Cancer-Gene-Expression
This repository includes a valuable source of datasets of clinical data of ovarian cancer gene expressions. 
Ovarian cancer is the most fatal gynecological malignancy among women. Making a reliable prediction of time to tumor recurrence would be a valuable contribution to post-surgery follow-up care. In this project we study three well-known data sets, known as TCGA, Tothill and Yoshihara, and compare three sparse regression methods, two of which (LASSO and Elastic Net) are well-known and the third (CLOT-combined l1 and l2 norm) is from our group at UT Dallas under the supervision of Professor Mathukumalli Vidyasagar. It is established that the three data sets are very different from each other. Therefore a two-stage predictor is built, whereby each test sample is first assigned to the most likely data set and then the corresponding predictor is used. The weighted concordance of each regression method is computed to compare the methods and select the best one. CLOT uses a biomarker panel of 103 genes and achieves a concordance index of 0.7829, which is higher than that achieved by the other two methods.

The detailed procedure can be found in our paper: https://link.springer.com/chapter/10.1007/978-3-319-59575-7_1
If you use the datasets, the code or the material in the paper, please cite the paper.




