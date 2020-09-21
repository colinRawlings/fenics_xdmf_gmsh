FROM quay.io/fenicsproject/stable:1.6.0

RUN apt-get update && apt-get install -y \
    gmsh \
    python-scipy \
    paraview

COPY ./requirements-dev.txt /home/fenics/requirements-dev.txt
RUN pip install -r /home/fenics/requirements-dev.txt

ENV LD_LIBRARY_PATH=/home/fenics/build/lib:$LD_LIBRARY_PATH
ENV PATH=/home/fenics/build/bin:/home/fenics/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
ENV PYTHONPATH=/home/fenics/build/lib/python2.7/site-packages:/home/fenics/legacy/pkgs:$PYTHONPATH
