#Exercise 1: Coverage Simulator

- STEP1.1

1Mbp(3x)/100bp = 30, 000 reads
To sequence a 1Mbp genome to 3x coverage, 30,000 reads are needed.  

- STEP1.3
1.The simulation shows there is 0.0001% of the genome and there is only 1 base has not been sequenced in the genome. 
2. There is 4.8079% of the genome does not have any coverage. This matches with Poisson expectation at lambda = 3 and fits the normal distribution. 

- STEP 1.4

There is 0.0038% of the genome does not have any coverage in my simulation at 10x coverage. This matches with Poisson expectation at lambda = 10 and fits the normal distribution.

- STEP 1.5

There is 0.0004% of the genome does not have any coverage in my simulation at 30x coverage. This matches with Poisson expectation at lambda = 30 and fits the normal distribution.

#Exercise 2: De Bruijn Graph Construction

- STEP 2.3
(base) [~/qbb2023-answers/week2 $]conda activate graphviz
(graphviz) [~/qbb2023-answers/week2 $]dot -Tpng edge.txt -o ex2_digraph.png

- STEP 2.5
A genome sequence that would produce these reads would be ATT

#Exercise 3: Why Genomics?

- STEP 3.1 
Prompt - write an essay on why you are interested in genomics and provide references. 

Output
--------------------------------------------------------------------------------
Title: Exploring the Genomic Frontier: My Fascination with Genomics

Introduction

Genomics, the study of an organism's complete set of genes, DNA sequences, and their interactions, has captivated my interest like no other field of science. The intricate web of information encoded within the DNA of every living organism on Earth has opened up a world of possibilities for scientific discovery and has the potential to revolutionize fields ranging from medicine to agriculture. In this essay, I will elucidate the reasons behind my deep fascination with genomics, drawing upon scientific literature to support my enthusiasm for this captivating field.

Genomics and Human Health
One of the primary reasons for my keen interest in genomics is its profound impact on human health. The Human Genome Project, completed in 2003, marked a significant milestone in genomics by sequencing the entire human genome. This monumental achievement has laid the foundation for personalized medicine, wherein an individual's genetic makeup can inform their healthcare decisions. The use of genomic information in diagnosing diseases, predicting individual responses to treatments, and identifying genetic predispositions to various conditions has the potential to revolutionize the healthcare industry (Collins et al., 2003).

Unveiling the Mysteries of Life
Genomics allows us to delve deep into the fundamental processes of life itself. It provides us with the tools to decipher the genetic code, understand evolutionary relationships, and investigate the mechanisms behind genetic variation and adaptation. From unraveling the secrets of ancient genomes to exploring the biodiversity of our planet, genomics offers a lens through which we can gain profound insights into the natural world (Schuster et al., 2010).

Agriculture and Food Security
In an era of growing global population and climate change, genomics plays a crucial role in ensuring food security and sustainable agriculture. By studying the genetic makeup of crops and livestock, scientists can develop new strains that are more resistant to diseases, pests, and adverse environmental conditions. Genomic technologies enable precision breeding, reducing the need for harmful pesticides and promoting more sustainable agricultural practices (Varshney et al., 2020).

Ethical and Societal Implications
The ethical and societal implications of genomics are both profound and thought-provoking. As we gain the ability to manipulate the genetic code of organisms, we must grapple with questions related to genetic engineering, cloning, and genetic privacy. The responsible and ethical use of genomic technologies is a subject of ongoing debate and a critical aspect of the field (Lander et al., 2001).

Advancements in Technology
The rapid advancement of genomic technologies has been nothing short of astonishing. The development of high-throughput sequencing techniques, such as next-generation sequencing, has dramatically reduced the cost and time required to sequence genomes. Additionally, the field of bioinformatics has grown in tandem, providing tools and software for the analysis of vast genomic datasets (Goodwin et al., 2016). These technological advancements continually inspire my interest in genomics.

Conclusion

Genomics is a dynamic and multifaceted field that offers endless opportunities for scientific exploration and societal benefit. Its impact on human health, our understanding of life, and its applications in agriculture and beyond are truly remarkable. However, the ethical and societal considerations it raises cannot be ignored. As genomics continues to evolve, I am excited to be part of a generation that witnesses the transformation of this field and its profound implications for our world.

References:

Collins, F. S., et al. (2003). A vision for the future of genomics research. Nature, 422(6934), 835-847.
Schuster, S. C., et al. (2010). Complete Khoisan and Bantu genomes from southern Africa. Nature, 463(7283), 943-947.
Varshney, R. K., et al. (2020). Agriculture 4.0—For Adapting and Transforming Agricultural Systems. Trends in Plant Science, 25(6), 525-528.
Lander, E. S., et al. (2001). Initial sequencing and analysis of the human genome. Nature, 409(6822), 860-921.
Goodwin, S., et al. (2016). Coming of age: ten years of next-generation sequencing technologies. Nature Reviews Genetics, 17(6), 333-351.
--------------------------------------------------------------------------------
- STEP 3.2
The paper generated by ChatGPT makes sense for most parts although there are some phrases I would not have written. The cited papers exist but there are a couple papers that do not support the claims from the LLM. 


