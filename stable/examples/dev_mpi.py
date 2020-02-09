#!python3

###############################################################
# Imports 
###############################################################

import dolfin as fn
from time import sleep

###############################################################
# Functions 
###############################################################

def log(comm: fn.MPI.comm_world, msg: str) -> None:
    print(f"log: rank {fn.MPI.rank(comm)}: {msg}")



###############################################################
# Main 
###############################################################

comm = fn.MPI.comm_world

print(
    f"mpi process: rank: {fn.MPI.rank(comm)}, "
    + f"size: {fn.MPI.size(comm)}")

log(comm, f"size: {fn.MPI.size(comm)}")

# try to use mpi barrier

# try to comm.bcast (a shared folder)