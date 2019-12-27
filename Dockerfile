FROM ceciledc/fenics_mixed_dimensional:latest

ARG vGMSH_VERSION=4.5.0

RUN pip install meshio
RUN pip install lxml
RUN pip install --no-binary=h5py h5py
RUN pip install jupyter
RUN pip install flake8 yapf jedi rope pytest

RUN mkdir /home/fenics/pkgs
COPY ./packages/fenics_utils /home/fenics/pkgs/fenics_utils
RUN pip install /home/fenics/pkgs/fenics_utils 

RUN apt-get update
RUN apt-get install -y gmsh libgfortran3
RUN mkdir /home/fenics/apps/
RUN curl -o /home/fenics/apps/gmsh.tgz http://gmsh.info/bin/Linux/gmsh-${vGMSH_VERSION}-Linux64-sdk.tgz
RUN cd /home/fenics/apps && tar -xvzf /home/fenics/apps/gmsh.tgz
ENV PYTHONPATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/lib:${PYTHONPATH}"
ENV PATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/bin:${PATH}"

CMD /usr/bin/bash