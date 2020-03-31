@page Initialization

We initialize M, A and B matrices for overlap alignment and non-overlap alignment.
* Aligning A<sub>i</sub> with B<sub>0</sub> without a gap. Not possible, hence, First column of M is initialized with -\f$\infty\f$.
* Aligning A<sub>0</sub> with B<sub>j</sub> without a gap. Not possible, hence, First row of M is initialized with -\f$\infty\f$.
* M<sub>0,0</sub> should have zero to begin the alignment using dynamic programming.

M | 0 | 1 | 2 | 3 | 4 | 5 |
--|--|--|--|--|--|---
0 | 0 | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ |
1 | -\f$\infty\f$ |  |  |  |  |  |
2 | -\f$\infty\f$ |  |  |  |  |  |
3 | -\f$\infty\f$ |  |  |  |  |  |

* Score of the best alignment between A<sub>i</sub> and B<sub>0</sub> that introduces a gap in A. Not possible, hence, First column of B is initialized with -\f$\infty\f$.
* Score of the best alignment between A<sub>0</sub> and B<sub>j</sub> that results a gap in B. Not possible, hence, First row of A is initialized with -\f$\infty\f$.

@section Init_section1 Overlap alignment

A | 0 | 1 | 2 | 3 | 4 | 5 |
--|--|--|--|--|--|---|
0 | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ |
1 | 0 |  |  |  |  | |
2 | 0 |  |  |  |  |  |
3 | 0 |  |  |  |  |  |


B | 0 | 1 | 2 | 3 | 4 | 5 |
--|--|--|--|--|--|---|
0 | -\f$\infty\f$ | 0 | 0 | 0 | 0 | 0 |
1 | -\f$\infty\f$ |  |  |  |  |  |
2 | -\f$\infty\f$ |  |  |  |  |  |
3 | -\f$\infty\f$ |  |  |  |  |  |

@section Init_section2 Non-overlap alignment

A | 0 | 1 | 2 | 3 | 4 | 5 |
--|--|--|--|--|--|---|
0 | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ | -\f$\infty\f$ |
1 | -22 |  |  |  |  | |
2 | -29 |  |  |  |  |  |
3 | -36 |  |  |  |  |  |


B | 0 | 1 | 2 | 3 | 4 | 5 |
--|--|--|--|--|--|---|
0 | -\f$\infty\f$ | -22 | -29 | -36 | -43 | -50 |
1 | -\f$\infty\f$ |  |  |  |  |  |
2 | -\f$\infty\f$ |  |  |  |  |  |
3 | -\f$\infty\f$ |  |  |  |  |  |
