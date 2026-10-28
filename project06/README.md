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
   ```r(i) = total distance of taxa i from all other taxa```
3. Calculate the new distance matrix (Q)
   ```Q(i,j) = d(i,j) - (r(i) + r(j))/(n-2)```
4. Use the new matrix (Q) to find the closest pair of taxa, i, j. Use the lowest distance and assign as the connecting node (u) for this pair.
5. Calculate branch length
  ```S(i,u) = d(i,j)/2 + (r(i)-r(j))/2(n-2)

     S(j,u) = d(i,j) - S(i,u)```
6. Calculate the new distances from (u) to the remaining taxa. The distance d(u,k) is
  ``` d(u,k) = [d(i,k) + d(j,k) - d(i,j)]/2 ```
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



# Struggles

*Github Repository*- 

*Biology Comprehension*-

*Connection of Functions*- 

# Personal Reflections
## Group Leader- Kailey 

- mention NP-complete projects as this one and the previous one were so. 

## Other member
Other members' reflections on the project

# Generative AI Appendix
Anthropic, Opus 5. 2026. Claude [Large language model]. https://claude.ai

Usage- Kailey 

Input: 

Output: 

Justification: 
