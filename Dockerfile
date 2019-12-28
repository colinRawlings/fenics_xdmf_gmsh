FROM ceciledc/fenics_mixed_dimensional:latest

ARG vGMSH_VERSION=4.5.0

RUN apt-get update

RUN apt-get install -y gmsh libgfortran3
RUN mkdir /home/fenics/apps/
RUN curl -o /home/fenics/apps/gmsh.tgz http://gmsh.info/bin/Linux/gmsh-${vGMSH_VERSION}-Linux64-sdk.tgz
RUN cd /home/fenics/apps && tar -xvzf /home/fenics/apps/gmsh.tgz
ENV PYTHONPATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/lib:${PYTHONPATH}"
ENV PATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/bin:${PATH}"

RUN pip install meshio
RUN pip install lxml
RUN pip install --no-binary=h5py h5py

RUN mkdir /home/fenics/pkgs
COPY ./requirements-dev.txt /home/fenics/pkgs/requirements-dev.txt
RUN pip install -r /home/fenics/pkgs/requirements-dev.txt