Multidevice Function and Options
--------------------------------
This code can run both on CPU (mostly for testing) and GPU.
For the gpu the following objects are relevant:

### `dev_info devinfo`
Contain the setup for the devices: how many ranks are there, which card to use, etc ...

### `pre_init_multidev1D`
Initialize the MPI library. Then the rank number and size is fetched, together with the processor name. 

The number of ranks is checked against the one given in the `geom_defines.txt`. If mismatched a panic is issued.