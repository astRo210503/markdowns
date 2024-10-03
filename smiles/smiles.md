
# SMILES:

overleaflink- 

## Gpt:
"SMILES" stands for Simplified Molecular Input Line Entry System. It's a notation system used to represent chemical structures in a way that can be easily interpreted by both humans and computers. Each atom in a molecule is represented by a symbol, and the bonds between atoms are represented by dashes, slashes, and colons.

On the other hand, "genome" typically refers to the complete set of DNA, including all of its genes. In a broader sense, it can refer to the complete genetic information of an organism.

GenAI can assist in various ways when it comes to generating and analyzing SMILES strings or genome sequences:

1. **Predicting Molecular Structures**: AI algorithms trained on vast datasets of chemical structures can predict the molecular structure based on a given SMILES string or other input.

2. **Drug Discovery**: GenAI can be utilized in drug discovery processes by generating and analyzing SMILES strings to identify potential drug candidates or optimize existing ones.

3. **Genetic Engineering**: In genetics, AI can help in generating synthetic DNA sequences for genetic engineering purposes, optimizing gene expression, or designing new organisms.

4. **Sequence Analysis**: AI algorithms can analyze genome sequences to identify patterns, mutations, or potential genetic disorders.

5. **Pharmacogenomics**: GenAI can analyze both genome sequences and drug response data to predict how individuals will respond to certain medications based on their genetic makeup.

Overall, genAI plays a crucial role in accelerating research and development in fields like pharmaceuticals, biotechnology, and personalized medicine by assisting scientists in generating, analyzing, and interpreting complex genetic and molecular data.

## Perplexity:

**What are SMILES (Chemistry)?**

SMILES (Simplified Molecular-Input Line-Entry System) is a notation system used to represent molecular structures in a compact and human-readable format. It is a line notation that uses a combination of letters and symbols to describe the atoms and bonds within a molecule. This system is widely used in chemistry and biochemistry to represent molecular structures in a concise and unambiguous manner, facilitating the exchange of information between researchers and the automation of various chemical processes.

**How can NLP Help in Generating SMILES?**

Natural Language Processing (NLP) can significantly contribute to the generation of SMILES by leveraging its capabilities in text analysis and machine learning. Here are some ways NLP can aid in generating SMILES:

1. **Text-based Molecular Structure Representation**: NLP can be used to analyze text-based representations of molecular structures, such as SMILES, and generate new molecular structures based on patterns and relationships identified in the text. This approach can be particularly useful for generating novel molecular structures that are not easily accessible through traditional methods[1][2][3][4].

2. **De Novo Molecule Generation**: NLP can be employed to generate new molecular structures by analyzing existing data points, extracting relevant information, and using this information to generate and optimize molecular structures based on desired properties. This approach can be particularly useful for hit identification in drug discovery, where potential drug candidates are identified from a large chemical space[2][4].

3. **Molecular Structure Optimization**: NLP can be used to optimize molecular structures by analyzing the properties of existing molecules and generating new structures that exhibit improved properties. This approach can be particularly useful for the design of novel drugs with enhanced efficacy and reduced side effects[2][4].

4. **Molecular Structure Prediction**: NLP can be employed to predict the properties of molecular structures based on their SMILES representation. This approach can be particularly useful for predicting the biological activity of novel molecules, which can significantly accelerate the drug discovery process[2][4].

In summary, NLP can play a crucial role in generating SMILES by leveraging its capabilities in text analysis and machine learning to analyze and generate molecular structures, optimize molecular properties, and predict molecular behavior. These applications have the potential to significantly accelerate and improve the efficiency of various chemical processes, including drug discovery[1][2][3][4].

Citations:
[1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9188322/
[2] https://www.nlpsummit.org/natural-language-processing-for-drug-discovery-the-state-of-practices-opportunities-and-challenges/
[3] https://www.linkedin.com/pulse/natural-language-processing-drug-discovery-aditya-gadiko
[4] https://developer.nvidia.com/blog/developing-nlp-applications-to-enhance-clinical-experiences-and-accelerate-drug-discovery/
[5] https://pubmed.ncbi.nlm.nih.gov/35702625/
[6] https://blog.gramener.com/nlp-use-cases-pharma-healthcare/
[7] https://openreview.net/pdf/VAVqG11WmSx0Wk76TAzp.pdf
[8] https://emerj.com/partner-content/nlp-drug-discovery-clinical/
[9] https://www.nature.com/articles/s41467-024-45102-8
[10] https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00642-3
[11] https://paperswithcode.com/task/text-based-de-novo-molecule-generation
[12] https://www.nature.com/articles/s41598-021-90259-7
[13] https://www.mdpi.com/2305-6304/12/2/153
[14] https://www.sciencedirect.com/science/article/abs/pii/S1359644620300465
[15] https://forums.fast.ai/t/create-language-model-for-chemical-structures/41398
[16] https://github.com/HUBioDataLab/SELFormer
[17] https://arxiv.org/abs/2304.04662
[18] https://towardsdatascience.com/basic-molecular-representation-for-machine-learning-b6be52e9ff76
[19] https://masterchemoinfo.u-strasbg.fr/en/blog/deep-learning-with-smiles

## Datasets:
### Hugging Face:
- maykcaldas/smiles-transformers (1 billion rows 35 GB) - https://huggingface.co/datasets/maykcaldas/smiles-transformers 
- HoangHa/belka-smiles-train-processed (98 million rows , 2.75 GB molecular smilles)- https://huggingface.co/datasets/HoangHa/belka-smiles-train-processed
- Models - https://huggingface.co/models?search=smiles

## Codes 
- ipynb
	- MolGAN (Generates new smiles) - https://github.com/deepchem/deepchem/blob/master/examples/tutorials/Generating_molecules_with_MolGAN.ipynb
- https://towardsdatascience.com/basic-molecular-representation-for-machine-learning-b6be52e9ff76
- https://github.com/deepchem/deepchem/blob/master/examples/tutorials/Introduction_to_GROVER.ipynb
---
[https://github.com/rsrinivas-repo/deepbind_molgen-](https://github.com/rsrinivas-repo/deepbind_molgen-) Repo for research project on generation of valid chemical SMILES using generative deeplearning algorithms.

  
  

[https://github.com/BadWrost/Smiles_generation](https://github.com/BadWrost/Smiles_generation) -Generate smiles base on max depth and specific atoms

  

[https://github.com/urchade/molgen](https://github.com/urchade/molgen) Molecule SMILES generation with GAN and Reinforcement Learning (Training Language GAN from scratch)

  
[https://github.com/shar032/Mol-Denoising-Autoencoder](https://github.com/shar032/Mol-Denoising-Autoencoder) to check invalid smiles


## References:
- https://www.perplexity.ai/search/what-are-smiles-oaltVmofS8y4cuWeNE4LHA
- https://github.com/seyonechithrananda/deepchem
- https://www.researchgate.net/publication/321374980_Deep_Reinforcement_Learning_for_De-Novo_Drug_Design
- https://github.com/EBjerrum/SMILES-enumeration/tree/master
- https://chat.mistral.ai/chat/2a6e0d8a-0f48-45a4-a542-c2aed4f1ff54
- https://chat.mistral.ai/chat/1fb42b85-2639-43d1-9070-8d6ec829b8fc
- https://dl.acm.org/doi/10.1145/3610533 imp
- https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2022.1046524/full imp
- https://github.com/CDDLeiden/DrugEx/tree/v3.2.0 imp
- https://github.com/LindeSchoenmaker/SMILES-corrector imp

## Python Libraries 
- deepchem - github - 
- rdkit
- pandas , numpy , matlplotlib

## Scratch files 
- Compound Generation with RDKit: Exploring Molecular Diversity for Given Elements and Total Number of Atoms - https://colab.research.google.com/drive/1FGaKuw06z_noH8nfyMa-H0NUHXPb0XlQ?usp=sharing
- https://github.com/urchade/molgen/blob/master/model.py

<!--stackedit_data:
eyJoaXN0b3J5IjpbNTA3MTU5NDA3LDE0NDgyODE5OTIsLTE2Nz
k4NDIyNTIsLTE2MTgyMzY0MjYsMTI3MTg2NzU5MSwzMTQ3NTU2
NjYsLTIxMzE4NDMxNjAsLTUxODQ2Mzg0MCwtMjY4NzkxNjgzLD
E4MzAzOTY1MDZdfQ==
-->
