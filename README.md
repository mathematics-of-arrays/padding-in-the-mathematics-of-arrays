# Padding in the Mathematics of Arrays 

This repository contains the code used in the experiments for the paper
[Padding in the Mathematics of Arrays](https://doi.org/10.1145/3460944.3464311),
as well as relevant references.

## Experiments

Experiments can be tuned and run for different hardware by adjusting
parameters within the Makefile in each subdirectory. The table below summarizes
which parameters are available and for which type of experiments:

| Parameter name | Description | Dimension lifting experiments | Padding experiments |
|----------------|-------------|-------------------------------|---------------------|
| *SIZE* | The length of each array axis (the experiments use 3-dimensional, cubic arrays) (default: 512) | yes | yes |
| *KSPLIT* | How to split the last axis of the arrays; for instance, a shape `(i, j, k)` whose last axis is split into *KSPLIT* slices has shape `(i, j,KSPLIT, k/KSPLIT)` (default: 4) | yes | no |
| *NTILES* | How to split all the axis of the arrays/tile the loops; a shape `(i, j, k)` that is split by into *NTILES* slices across each axis has shape `(NTILES, i/NTILES, NTILES, j/NTILES, NTILES, k/NTILES)` (default: 4) | yes | no |
| *THREADS* | How many threads to use for computation (default: 4) | yes | no |


### Baseline

The original snippet we use for our baseline experiments is implemented in
`common/original.c`.

### Dimension Lifting Experiments

Building the dimension lifting experiments produces four binaries, all
different versions of the PDE solver snippet studied in the paper
[Finite difference methods fengshui: alignment through a mathematics of
arrays](http://doi.org/10.1145/3315454.3329954):

- `bin/S.bin`, which is built from `common/original.c` and serves as a baseline
  for the experiments;
- `bin/MDL.bin`, in which the arrays are lifted into *THREADS* slices along
  their first axis; program execution is then parallelized over *THREADS*
  threads along the resulting first axis;
- `bin/MDLSL.bin`, in which the arrays are lifted into *KSPLIT* slices along
  their third axis; program execution is then parallelized over *THREADS*
  threads along the resulting third axis;
- `bin/MDLTM.bin`, in which the array is lifted into *NTILES* slices along each
  dimension. The loop structure is tiled; the outer loops are collapsed, and
  program execution is then parallelized over *THREADS* threads.

### Padding experiments

Building the padding experiments produces three binaries, all different versions
of the PDE solver snippet studied in the paper [Finite difference methods
fengshui: alignment through a mathematics of
arrays](http://doi.org/10.1145/3315454.3329954):

- `bin/nopad.bin`, which is built from `common/original.c` and serves as a
  baseline for the experiments;
- `bin/padj.bin`, in which the arrays are *circularly* padded along their
  second axis once on each side at. Padding is replenished at the end of each
  iteration along the second axis;
- `bin/padk.bin`, in which the arrays are *circularly* padded along their
  third axis once on each side at. Padding is replenished at the end of each
  iteration along the third axis.

## References

We have started implementing the Mathematics of Arrays formalism using Coq. The
(draft) code is available in the
[mathematics-of-arrays/moa-formalization](https://github.com/mathematics-of-arrays/moa-formalization)
GitHub repository.
