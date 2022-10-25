Control Flow
------------

This file broadly document how the program run, following the control flow.

### `main` in `src\main.c`:

#### Reading setups

The computer time is fetched, and stored in the global structure `mc_params`.

A message is logged on which gauge action is used (see `GAUGE_ACTION` in [constants.md])

If the compilation target a `MULTIDEVICE` machine, the initialization procedures `pre_init_multidev1D` e `gdbhook` are runned.

After some logging (done only by the parent process) the main input procedure `set_global_vars_and_fermions_from_input_file` is called. 