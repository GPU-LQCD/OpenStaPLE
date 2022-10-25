Support for Gnu DeBugger
------------------------

Some support for GDB is given by the function `void gdbhook()`:
if the enviroment variable `GDBHOOK` evaluate to non-zero a volatile flag is set to 0 and waited for.
The process are then sincronized with `MPI_Barrier` and the execution resume