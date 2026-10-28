# Introduction

--- 

## Neighbor-joining: Building Unrooted Phylogenetic Trees

Neighbor-joining (NJ) is a greedy agglomerative clustering algorithm used to construct unrooted phylogenetic trees from a distance matrix derived from the results of an alignment score matrix. It is used to compute genetic differences between species or variants and analyze where they may have diverged at the molecular level. Unlike UPGMA (Unweighted Pair Group Method with Arithmetic Mean), neighbor-joining:


1. Produces an unrooted tree
2. Removes the assumption that evolution occurs at a constant rate
3. Uses an additive metric

Each of these represents a unique characteristic of the neighbor-joining method.

# Algorithm

Within this project, the Neighbor-joining algorithm is followed by the local alignment algorithm, Smith-Waterman. The NJ algorithm steps root from the guide provided in the Jupyter Notebook as we constructed our function around this: 

**Initialization**:
- Let $n$ be the number of sequences
- Assign each sequence $i$ to its own leaf node
- Initialize tree $T$ with leaf nodes

**Iteration**:
While $n > 2$:

Pre-step: Calculate the total distance (TD) of each row in the matrix 

1. Calculate Q-matrix (D*) where:
   $Q(i,j) = (n-2)d(i,j) - \sum_{k=1}^n d(i,k) - \sum_{k=1}^n d(j,k)$

2. Find pair $(i,j)$ with minimum $Q(i,j)$

3. Calculate branch lengths:
   $d_i = \frac{d(i,j) + (r_i - r_j)/(n-2)}{2}$
   $d_j = d(i,j) - d_i$
   where $r_i = \sum_{k=1}^n d(i,k)$

4. Create new node $k$
   - Add branches from $k$ to $i$ and $j$ with lengths $d_i$ and $d_j$

5. Update distances to remaining nodes $x$:
   $d(k,x) = \frac{d(i,x) + d(j,x) - d(i,j)}{2}$

6. Remove nodes $i$ and $j$
7. Add node $k$ to active nodes
8. $n = n - 1$

**Termination**:
When $n = 2$ with remaining nodes $i$ and $j$:
- Add final branch between $i$ and $j$ with length $d(i,j)$
- Return unrooted tree $T$


# Pseudocode
Put pseudocode in this box:

```
Some pseudocode here
```

## Biological Considerations and Context

To expand on neighbor-joining, several variants have been developed since the original algorithm:

1. Live neighbor-joining
2. Extended neighbor-joining
3. Multifurcating neighbor-joining

Each variant builds upon the original algorithm and introduces additional parameters to accommodate different types of data.

One area where our group experienced the most confusion was the conversion of sequences into genetic distances. Based on the provided HTML, the process begins by calculating a p-distance. However, p-distance can underestimate evolutionary change as sequence divergence increases. Therefore, an important consideration is choosing between the Jukes-Cantor and Kimura 2-parameter models.

* **Jukes-Cantor**
    * Assumption: All four nucleotides occur with equal frequencies, and all nucleotide substitutions occur at the same rate.

* **Kimura 2-Parameter (K2P)**
    * Assumption: All four nucleotides occur with equal frequencies, but transitions occur more frequently than transversions. The model uses two substitution parameters: one for transitions and one for transversions.


Both models were challenging for us to fully understand. With additional time, we would be better able to determine how each model operates and which would be most appropriate for this dataset.


The dataset consists of HIV sequences. An important biological consideration is the complexity of viral evolution. High mutation rates and recombination events can lead to inaccurate branching patterns and should be considered when constructing a phylogenetic tree.

### Strengths

- Computationally efficient
- Low computational cost, O(N³)
- Well suited for smaller datasets, such as the viral dataset used in this project

### Limitations

A few limitations of neighbor-joining include:

- Can assign negative branch lengths in some cases
- Results can be difficult to reproduce depending on model parameters and implementation details
- Like other clustering methods, accuracy is highly dependent on the quality of the distance matrix

### Tools to Compare

We unfortunately did not have the time to generate comparison between a known Bioinformatics tool due to time constraints. This is the tool we would have hoped to utilize that could provide a comparison of our results:

https://www.hiv.lanl.gov/components/sequence/HIV/treemaker/treemaker.html

This one is interesting because it uses R rather than Python to emulate the neighbor joining algorithm and purposefuly is marked as the HIV sequence database. 

# Successes

* *Constructing a full pipeline*- Built the full pipeline, reading the FASTA file, normalized Smith–Waterman scoring, a symmetric 20 × 20 distance matrix, and neighbor joining that outputs a Newick tree.
* *Distance Transformation*: Chose and justified a distance transformation. We normalized scores by the geometric mean of the self-scores and used distance = 1 − similarity. This gives 0 for identical sequences, is symmetric, and is never negative.
* *Validation*- We successfully verified our code with known answers. The neighbor-joining code reproduced the correct trees for a 4-taxon example, the walkthrough’s example, and the additive HIV matrix.
* *Review*- We caught problems in the provided materials: the docstring expects 19 sequences but the file has 20, the headers needed parsing to get usable names, and two worked values in the walkthrough were miscalculated.
*Full Deep-dive into Numpy*- We learned new tools from the Numpy module such as NumPy’s delete, and other functions that we can use to facilitate coding with using arrays/ matrices. Using the Numpy user guide provided a useful resource with implementation and viewing examples of how it can be used. 

# Struggles

* *Content Overload*- Reconciling two sources was confusing for us at first. The notebook and the walkthrough use different letters for the same things (Q vs. D*, rᵢ vs. TD, k vs. u) and stop at different points. Switching between them caused confusion until we picked one and mapped the other onto it.

* *Series of Formulas*- It was also a bit of a hustle turning the formulas into code. We split the branch-length formula across two lines and accidentally counted the distance twice. We also left out a / 2. Neither caused an error message. They just gave wrong numbers. These were corrected.

* *Algorithm Comprehension*- Understanding the algorithm proved challenging due to the number of steps and formulas involved. The HTML resource was helpful for organizing the workflow and providing a simple example. The most time-consuming aspect was developing a thorough understanding of the concepts before implementing the code.

* *Building the Distance Matrix from Alignment Scores*- This stage was both successful and challenging. We encountered conflicting information when determining how to compute distances and construct the matrix. Both Jukes-Cantor and K2P have distinct strengths and limitations. Given more time, we would like to revisit this topic to better evaluate the criteria for selecting one model over the other.

# Personal Reflections
## Group Leader- Kailey 

- mention NP-complete projects as this one and the previous one were so.
- mention using R vs Python

## Other member
Other members' reflections on the project

Dianah:
This project was a bit challenging, and it took us a while to figure out. After our first virtual meeting, we each went back to review the slides, notebook and walkthrough, but even at our second meeting we were still struggling to see how the pieces fit together.It taught me that understanding an algorithm and implementing it are two different skills. The team discussions helped a lot. I was able to follow neighbor joining on paper well but writing it correctly was confusing. Python also accepted some small bugs without complaint which quietly produced wrong answers. The notebook, the walkthrough and the data file each had small inconsistencies, and checking them against each other made me more confident in my own results.

# References 

Segura-Alabart, N., Serratosa, F., & Fernández, A. (2024). A practical study of the proportion of non-unique neighbor-joining trees of microsatellite markers. Computational and Structural Biotechnology Reports, 1, Article 100022. https://doi.org/10.1016/j.csbr.2024.100022


# Generative AI Appendix
Anthropic, Opus 5. 2026. Claude [Large language model]. https://claude.ai

Usage- Kailey 

Input: How is the Smith-Waterman algorithm connect to the NJ algorithm? How to approach building a distance matrix?

Justification: I used Claude to provide assistant in detailing how the concepts are integrated together. It provided me resources to look at how to compute the Distance Matrix and develop pseudocode for it. 

Dianah
I used Claude to learn the use of np.delete, np.vstack and np.column_stack to remove the joined rows and columns and add the new node’s row and column. I also used it to understand the concepts of the project.
