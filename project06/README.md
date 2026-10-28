# Introduction

--- 

## Neighbor Joining: Building Unrooted Phylogenetic Trees

Neighbor Joining is a greedy heuristic clustering algorithm used to construct unrooted phylogenetic trees from a distance matrix crafted by the results of an alignment score matrix. This is used to compute genetic differences between species or variations and provide analysis of where they might have diverged at the molecular level. Unlike the UPGMA (Unweighted Pair Group Method with Arithmetic Mean), neighbor-joining 

1. Produces an unrooted tree
2. Removes the assumption that evolution occurs at a constant rate
3. Uses an additive metric

Each provide as a unique aspect to neighbor joining. 

# Algorithm (Note, I will match the variables to our code)

Within this project, the Neighbor Joining algorithm is followed by the local alignment algorithm, Smith-Waterman. The algorithm steps proposed here are what happens after given the alignment scores from Smith-Waterman. 

1. Convert the alignment score matrix into a distance-based matrix
2. Calculate the net divergence for each taxa (across the rows)
   r(i) = total distance of taxa i from all other taxa
3. Calculate the new distance matrix (Q)
   Q(i,j) = d(i,j) - (r(i) + r(j))/(n-2)
4. Use the new matrix (Q) to find the closest pair of taxa, i, j. Use the lowest distance and assign as the connecting node (u) for this pair.
5. Calculate branch length
   S(i,u) = d(i,j)/2 + (r(i)-r(j))/2(n-2)

   S(j,u) = d(i,j) - S(i,u)
6. Calculate the new distances from (u) to the remaining taxa. The distance d(u,k) is
   d(u,k) = [d(i,k) + d(j,k) - d(i,j)]/2
7. Make steps 2-6 recursive until n > 2
8. 

# Pseudocode
Put pseudocode in this box:

```
Some pseudocode here
```

## Biological Considerations and Context

### Strengths

- Computationally fast

### Limitations

A few limitations to this project are

- Assigns negative length to some branches often.
-

### Tools to Compare

We unfortunately did not have the time to generate comparison between a known Bioinformatics tool due to time constraints. This is the tool we would have hoped to utilize that could provide a comparison of our results:

https://www.hiv.lanl.gov/components/sequence/HIV/treemaker/treemaker.html

This one is interesting because it uses R rather than Python to emulate the neighbor joining algorithm and purposefuly is marked as the HIV sequence database. 

# Successes
* Built the full pipeline, reading the FASTA file, normalized Smith–Waterman scoring, a symmetric 20 × 20 distance matrix, and neighbor joining that outputs a Newick tree.
* Chose and justified a distance transformation. We normalized scores by the geometric mean of the self-scores and used distance = 1 − similarity. This gives 0 for identical sequences, is symmetric, and is never negative.
* We successfully verified our code with known answers. The neighbor-joining code reproduced the correct trees for a 4-taxon example, the walkthrough’s example, and the additive coronavirus matrix.
* We caught problems in the provided materials: the docstring expects 19 sequences but the file has 20, the headers needed parsing to get usable names, and two worked values in the walkthrough were miscalculated.
* We learned new tools, NumPy’s delete, vstack and column_stack, building Newick strings, and debugging from tracebacks.


# Struggles
* Reconciling two sources was confusing for us at first. The notebook and the walkthrough use different letters for the same things (Q vs. D*, rᵢ vs. TD, k vs. u) and stop at different points. Switching between them caused confusion until we picked one and mapped the other onto it.
* It was also a bit of a hustle turning the formulas into code. We split the branch-length formula across two lines and accidentally counted the distance twice. We also left out a / 2. Neither caused an error message. They just gave wrong numbers. These were corrected.

*Github Repository*- 

*Biology Comprehension*-

*Connection of Functions*- 

# Personal Reflections
## Group Leader- Kailey 

- mention NP-complete projects as this one and the previous one were so. 

## Other member
Other members' reflections on the project

Dianah:
This project was a bit challenging, and it took us a while to figure out. After our first virtual meeting, we each went back to review the slides, notebook and walkthrough, but even at our second meeting we were still struggling to see how the pieces fit together.It taught me that understanding an algorithm and implementing it are two different skills. The team discussions helped a lot. I was able to follow neighbor joining on paper well but writing it correctly was confusing. Python also accepted some small bugs without complaint which quietly produced wrong answers. The notebook, the walkthrough and the data file each had small inconsistencies, and checking them against each other made me more confident in my own results.

# Generative AI Appendix
Anthropic, Opus 5. 2026. Claude [Large language model]. https://claude.ai

Usage- Kailey 

Input: 

Output: 

Justification: 

Dianah
I used Claude to learn the use of np.delete, np.vstack and np.column_stack to remove the joined rows and columns and add the new node’s row and column. I also used it to understand the concepts of the project.
