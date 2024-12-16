

## pole class

pole class has the following attributes:

- position
- weights
- normal vector
- updater function (to update the intensity, that is the potential, of the pole)

## Buckets class

Buckets class stores specified objects as `Buckets<T>`, and generates tree structure until the number of objects in a bucket is less than or equal to the specified number of objects per bucket.

The step to generate the tree structure should be as follows:

1. add objects to the bucket
2. set the maximum level of the tree using `setLevel`
3. generate the tree structure using `generateTree` while specifying the condition to stop the generation of the tree structure


# Fast Multipole Method

The Fast Multipole Method (FMM) is an algorithm for the efficient calculation of the integration of the pole/potential using the tree structure, the multipole expansion, shifting expansion, and the local expansion. Since FMM calculates integration/summation, such as BIE and does not make the coefficient matrix, solver for the simultaneous linear equations should be iterative methods. GMRES is commonly used for the solver with FMM.

| First steps | GRMES iterative step | description | | |
| --- | --- | --- | --- | --- |
| 1 | | add poles to the root bucket | | |
| 2 | | generate the tree structure from the root bucket | | |
| 3 (before M2M) | | expansion of the poles | | |
| 4 | 1 | **update the intensity of the poles** | | |
| 5 | 2 | Multipole to Multipole (M2M): shift the multipole expansion at each center, from the deeper level to the upper level | about 8 🪣 -> 1 parent 🪣 | use pre-computed SPH |
| 6 | 3 |  Multipole to Local (M2L)| every 🪣 -> (only same level) -> many local 🪣 | use pre-computed SPH |
| 7 | 4 | Local to Local (L2L) | 1 🪣 -> about 8 children 🪣 | use pre-computed SPH |
| 8 | 5 | Add direct integration for the near field and the integration using the local expansion for the far field | | |

Many part of process are dependent on relative position of the poles and the buckets. Therefore, many part of the first steps are saved and reused in the following iterative steps. Remaining part for iterative steps are the update of the intensity of the poles, and simple incrementatation in four-fold for-loops. However, the number of incrementation is not negligible, and the direct integration for the near field also takes time. FMM is surely faster than the direct summation when the number of poles is more than about 10000, but the calculation time is already long when the number of poles is about 10000.