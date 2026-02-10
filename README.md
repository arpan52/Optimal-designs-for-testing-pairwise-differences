# Optimal designs for testing pairwise differences: a graph based game theoretic approach

# Description
This repository contains all the original **R** codes used for the numerical studies in the paper entitled **Optimal designs for testing pairwise differences: a graph based game theoretic approach** authored by **_Arpan Singh, Satya Prakash Singh and Ori davidov_**. For a through reading, refer to the [Paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12757).

# Instructions
 To reproduce the results, use the original **R** codes that are available in the **Codes** branch. 

1) The following **R** packages are required to use the codes,\
      a) **nloptr**\
      b) **mvtnorm**
      
2) Names of the main files and the **R** codes are self-explainatory. For example the file **Table 2 Supplement (UIT Path)** contains all the **R** codes used to produce the results given in the **Table 2** in the **Supplementary material**, in which **Path** comparison is done for the **UITs**. 

3) Similarly, the **R** code **UIT_tree_K3.r**  produce the max--min design for the **Tree** graph for **K=3 groups** when **UIT** test is used. Note that **Tree Graph** has been renamed as **Star Graph** in the latest version of the paper. Thus, analogously, all the codes that contain the word **Tree** should be understood as the code of a **Star** graph.

4) To reproduce the results of Table 3, use the new **R** codes in the file **Simulations empirical power Table 3** instead of the old version **Table 3 (Simulations real data)**.


