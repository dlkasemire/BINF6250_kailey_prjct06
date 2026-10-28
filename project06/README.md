# Introduction

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
def read_fasta(filename: str) -> Dict[str, str]:
initialize empty dict seqs, current_id = None
        open file, for each line:
            strip whitespace, skip if empty
            if line starts with ">":
                current_id = last token of header - after ">"
                seqs[current_id] = ""
            else if current_id exists:
                append line to seqs[current_id]  (handles multi-line sequences)
        return seqs

def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap, best_score):
if i == 0 or j == 0:
            return 0, best_score  (first row/column are boundary cells)
        diag_score = matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
        up_score   = matrix[i-1][j] + gap
        left_score = matrix[i][j-1] + gap
        score = max(0, diag_score, up_score, left_score)  (0 = start new local alignment)
        best_score = max(best_score, score)
        return score, best_score

def smith_waterman(
    sequence1: str,
    sequence2: str,
    match_score: float = 2.0,
    mismatch_penalty: float = -1.0,
    gap_penalty: float = -1.0
) -> float:
  initialize (len1+1) x (len2+1) score matrix of 0s, best_score = 0
    for each cell (i, j):
        score_matrix[i][j], best_score = cal_score(...)
    self1 = match_score * len(sequence1)  (score of seq1 aligned to itself)
    self2 = match_score * len(sequence2)
    similarity = best_score / sqrt(self1 * self2)  (geometric mean normalizes to 0-1)
    return similarity  (1.0 = identical; distance conversion done in build_distance_matrix)

def build_distance_matrix(sequences: Dict[str, str]) -> Tuple[np.ndarray, List[str]]:
seq_ids = list of sequence IDs  (fixes row/column order)
        n = number of sequences
        initialize n x n matrix of 0s  (diagonal stays 0)
        for each pair i < j:
            similarity = smith_waterman(seq_i, seq_j)
            distance = 1 - similarity
            matrix[i, j] = matrix[j, i] = distance  (symmetric)
        return matrix, seq_ids

def neighbor_joining(
    distance_matrix: np.ndarray,
    labels: List[str]
) -> str:
copy distance matrix, nodes = copy of labels, n = number of nodes
        if n == 0: return ";"    if n == 1: return "label;"
        while n > 2:
            1. row_totals TD(i) = sum of row i
            2. Q(i,j) = (n-2) * D(i,j) - TD(i) - TD(j)
               set diagonal of Q to infinity  (a node can't join itself)
            3. (i, j) = position of minimum Q  (i < j since Q is symmetric)
            4. branch lengths:
               di = D(i,j)/2 + (TD(i) - TD(j)) / (2(n-2))
               dj = D(i,j) - di
               clamp negative branch lengths to 0
            5. new node label = "(node_i:di,node_j:dj)"  (Newick string builds itself)
            6. distances to new node: D(u,k) = (D(i,k) + D(j,k) - D(i,j)) / 2
            7. overwrite row/column i with new distances, set D(u,u) = 0
               delete row and column j
               nodes[i] = new label, delete nodes[j]
            8. n = n - 1
        when 2 nodes remain:
            split final distance d in half between them
            return "(node_0:d/2,node_1:d/2);"

```

## Biological Considerations and Context

### Expansions of Neighbor-joining

To expand on neighbor-joining, several variants have been developed since the original algorithm:

1. Live neighbor-joining
2. Extended neighbor-joining
3. Multifurcating neighbor-joining

Each variant builds upon the original algorithm and introduces additional parameters to accommodate different types of data.

### Conversion of Sequences into Genetic Distances

One area where our group experienced the most confusion was the conversion of sequences into genetic distances. Based on the provided HTML, the process begins by calculating a p-distance. However, p-distance can underestimate evolutionary change as sequence divergence increases. Therefore, an important consideration is choosing between the Jukes-Cantor and Kimura 2-parameter models.

* **Jukes-Cantor**
    * Assumption: All four nucleotides occur with equal frequencies, and all nucleotide substitutions occur at the same rate.

* **Kimura 2-Parameter (K2P)**
    * Assumption: All four nucleotides occur with equal frequencies, but transitions occur more frequently than transversions. The model uses two substitution parameters: one for transitions and one for transversions.


Both models were challenging for us to fully understand. With additional time, we would be better able to determine how each model operates and which would be most appropriate for this dataset. As a note, we opted for a normalization method where we computed 1 - Smith-Waterman similarity. 

### The Data itself 

The dataset consists of HIV sequences. An important biological consideration is the complexity of viral evolution. High mutation rates and recombination events can lead to inaccurate branching patterns and should be considered when constructing a phylogenetic tree.

## Strengths

- Computationally efficient
- Low computational cost, O(N³)
- Well suited for smaller datasets, such as the viral dataset used in this project

## Limitations

A few limitations of neighbor-joining include:

- Can assign negative branch lengths in some cases
- Results can be difficult to reproduce depending on model parameters and implementation details
- Like other clustering methods, accuracy is highly dependent on the quality of the distance matrix

## Tools to Compare

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

* *Building the Distance Matrix from Alignment Scores*- This stage was both successful and challenging. We encountered conflicting information when determining how to compute distances and construct the matrix. Both Jukes-Cantor and K2P have distinct strengths and limitations, yet we opted for a 1 - Smith-Waterman normalization approach. Given more time, we would like to revisit this topic to better evaluate the criteria for selecting one model over the other.

# Personal Reflections
## Group Leader- Kailey 

This project was not an easy one to fathom. The algorithm was truly a difficult one to comprehend with the various formulas and steps. It's different learning about phylogenetic trees from biology and noting that the different branches are different lengths to represent divergence in evolutionary relationships. This project said "how do we build those phylogenetic trees." This was never a consideration I had until this course. 

Completing the HTML provided to be an extraordinary resource as it provided a walk through of the entirety of the algorithm. I was able to do the example provided by hand to gain a better understanding of the material. Reading scientific papers and additional resources also honed down on how these algorithms are used. What was interesting is that a lot more papers mentioned UPGMA rather than neighbor-joining so I would like to continue my search into why that is the preference. 

We were able to meet a couple times throughout this process, and for our initial meeting we had spent the majority of the time processing the material. The second meeting we had code about halfway finished, but the discussion persisted about confusion we had about the concept and formulas. Peer collaboration was a useful asset here. I would also like to thank my group for also being accommodating of the West Coast schedule because there were a few times where we would meet and talk late at night! 

## Dianah
This project was challenging, and it took us a while to figure out. We went through the project together in our first meeting, but even after this, we still needed to understand the concepts some more. We each went back to review the slides, notebook and walkthrough, which made our second meeting much better although we were still struggling to see how all the pieces fit together. The team discussions helped a lot. I was able to follow neighbor joining on paper well but writing it correctly was confusing. The notebook, the walkthrough and the data file each had small inconsistencies, and python also accepted some small bugs without complaint which quietly produced wrong answers. All this taught me that understanding an algorithm and implementing it are two different skills.

## Yulia
This project was definitely a challenge as there were many working parts. First we had to get a good grasp of how neighbor joining works. Our group met on Teams call and spent the first meeting just trying to understand the project, which took so long that we only started working on it actually in our second meeting. In that first meeting, we went through the NJ html activity together to understand the formulas and repetition while n>2, and the termination steps at n=2 with the distance matrices. In the second meeting, we discussed code ideas for implementing Smith-Waterman to build the initial distance matrix, and preliminary code for the NJ function. Using TODO blocks (Kailey's suggestion) was very helpful as there were definite gaps in our code that we knew we had to fill, but we did not know how to phrase it yet. Overall, I feel I have a good understanding of the concepts of Neighbor Joining now, although I still remain confused on if our method of calculating distances from similarities was the most optimal. 


# References 

Segura-Alabart, N., Serratosa, F., & Fernández, A. (2024). A practical study of the proportion of non-unique neighbor-joining trees of microsatellite markers. Computational and Structural Biotechnology Reports, 1, Article 100022. https://doi.org/10.1016/j.csbr.2024.100022


# Generative AI Appendix
Anthropic, Opus 5. 2026. Claude [Large language model]. https://claude.ai

Usage- Kailey 

Input: How is the Smith-Waterman algorithm connect to the NJ algorithm? How to approach building a distance matrix?

Justification: I used Claude to provide assistant in detailing how the concepts are integrated together. It provided me resources to look at how to compute the Distance Matrix and develop pseudocode for it. 

--- 
Dianah
Claude was very helpful in explaining the concepts of the assignment and working through parts of the code that were new to me.

-------------------
Usage: Yulia

Input: What are good methods for converting Smith-waterman similarity score to distances?

Justification: Although we had good understanding of Smith-Waterman for similarities and Neighbor Joining, we were confused on how to make the distance matrix. Using the Claude output, we decided to use 1 − normalized Smith-Waterman similarity = distance.

-------------------
Usage: Yulia

Input: Update the pseudocode to match the updated code

Justification: After working through the code together, it deviated from our original pseudocode so this helped make sure the pseudocode was up to date.
