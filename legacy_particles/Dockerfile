FROM quay.io/fenicsproject/stable:1.6.0

RUN apt-get update && apt-get install -y gmsh 
RUN sudo apt-get install python-scipy

# ENV PYTHONPATH="/home/fenics/pkgs:/home/fenics/build/lib/python2.7/site-packages"

COPY ./requirements-dev.txt /home/fenics/requirements-dev.txt
RUN pip install -r /home/fenics/requirements-dev.txt

# TODO: cross platform strategy !!!
ENV DISPLAY=192.168.0.22:0

ENV LD_LIBRARY_PATH=/home/fenics/build/lib:$LD_LIBRARY_PATH
ENV PATH=/home/fenics/build/bin:/home/fenics/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
ENV PYTHONPATH=/home/fenics/build/lib/python2.7/site-packages:$PYTHONPATH
