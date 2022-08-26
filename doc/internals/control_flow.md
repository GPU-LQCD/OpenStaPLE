Control Flow
------------

This file broadly document how the program run, following the control flow.

### `main` in `src\main.c`:

#### Reading setups

The computer time is fetched, and stored in the global structure `mc_params`.

A message is logged on which gauge action is used (see `GAUGE_ACTION` in [constants.md])