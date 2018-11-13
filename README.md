# Multi Locus Analysis

A collection of utilities for analyzing multiple loci (diffusing).

### Philosophy

The routines tend to be centered around Pandas, attempting to be general by
allowing you to specify which columns to groupby and which columns to perform
each operation on freely.

Whenever Pandas leaks memory, or sane groupby/apply semantics would lead to
impossibly large intermediate dataframes, routines will instead dump these
intermediate results to disk, where they can be read in again later.

Names and conventions reflect the fact that most calculations that you could
want to do to analyze trajectories take the form of calculations of either
distributions or moments of some discrete derivative of the motion. (e.g. MSD is
1st moment, 1st derivative; Velocity Cross-Correlation is 2nd moment, first
derivative; displacement or distance distributions are 1st or 0th derivative,
respectively...).


