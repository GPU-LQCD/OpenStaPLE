#include "../Mpi/multidev.h"
#include <unistd.h>
#include <mpi.h>
#include "debugger_hook.h"


void gdbhook()
{
    volatile int flag = 0;
    if(0==devinfo.myrank && getenv("GDBHOOK")){
        printf("MPI%02d: Thank you for using the GDB HOOK!\n",devinfo.myrank);
        printf("         Waiting for user intervention.Please attach to process %d,\n",
                getpid());
        printf("         [hint: gdb -p %d ] \n",  getpid());
        printf("         and set the value of 'flag' to 1.\n");
        printf("         You may have to use 'finish' a couple of times to go down the call stack.\n");
        printf("         HAPPY HUNTING AND GOOD LUCK!!!!\n");
        while(1 != flag)
            sleep(1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
