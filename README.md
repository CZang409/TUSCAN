# TUSCAN

`TUSCAN` is an R package that aims to perform tumor segmentation and classification analysis by constructing a spatial copy number profile. TUSCAN takes both the gene information from spatial transcriptomics data and the hematoxylin and eosin (H&E) staining image as input to improve the accuracy of tumor region identification. Below is the overall workflow of TUSCAN, which consists of three main steps. The first step is to find a subset of normal spots. TUSCAN performs a preliminary
clustering of all spots based on the gene information dimensionally reduced by an autoencoder, along with RGB values extracted from the histology image. Then, it will automatically pick a cluster with the highest confidence of representing normal cells. The second step is to infer the copy number profile of all spots applying the selected cluster as the reference baseline. The third step is to segment all spots into tumor regions and normal regions via a consensus clustering with the copy number profile as the input.

![Untitled (3)](https://github.com/CZang409/TUSCAN/assets/166551317/cdccc0c6-6feb-47ce-9782-f36a044eae2e)
