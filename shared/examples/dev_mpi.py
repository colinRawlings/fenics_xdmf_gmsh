#!python3

###############################################################
# Imports 
###############################################################

import dolfin as fn
from time import sleep
from datetime import datetime

###############################################################
# Functions 
###############################################################

def log(comm: fn.MPI.comm_world, msg: str) -> None:
    print(f"{datetime.now()}: log: rank {fn.MPI.rank(comm)}: {msg}")

###############################################################
# Main 
###############################################################

comm = fn.MPI.comm_world

print(
    f"mpi process: rank: {fn.MPI.rank(comm)}, "
    + f"size: {fn.MPI.size(comm)}")

log(comm, f"size: {fn.MPI.size(comm)}")

# test barrier

if fn.MPI.rank(comm) == 0:
    sleep_time = 0
else:
    sleep_time = 1

log(comm, f"starting to sleep")
sleep(sleep_time)

log(comm, f"done sleeping, waiting at barrier")
comm.Barrier()
log(comm, f"past barrier")

# test communication

if fn.MPI.rank(comm) == 0:
    tmp_filepath = "fil142141.log"
    log(comm, f"init tmp_filepath: {tmp_filepath}")
else:
    tmp_filepath = None
    log(comm, f"init tmp_filepath: {tmp_filepath}")


tmp_filepath = comm.bcast(tmp_filepath, root=0)

log(comm, f"tmp_filepath: {tmp_filepath}")


# try to use mpi barrier
#
# try to comm.bcast (a shared folder)