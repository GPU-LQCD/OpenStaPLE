Preprocessor Constants
---------

This file collect preprocessor constant, what they seemily do and where are used.

### Physical
Values that are part of the physical simulation.

#### `GAUGE_ACTION`
Check which gauge is used:
- 1      : Three-level Symanzik improved gauge action
- default: Wilson gauge action

### Target specification
Values that specify for what target the program is compiling for.

#### `MULTIDEVICE`
If the code is parallelized on a GPU or runned on a simple CPU.
See [multidevice.md]

#### `NRANKS_D3`
Number of ranks to use. `D3` refers to the *salamino* model, meaning the lattice should be divided in `NRANKS_D3` equal slices.

### Others
Values that are not part of any category.

#### `COMMIT_HASH`
The current git head hash. Collected by `configure.ac`.