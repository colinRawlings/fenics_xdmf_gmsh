FROM ceciledc/fenics_mixed_dimensional:latest 
# one container older: v2019.1

ARG vGMSH_VERSION=4.5.0
ARG vPYRIGHT_VERSION=1.1.19

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

RUN curl -sL https://deb.nodesource.com/setup_12.x | sudo -E bash -
RUN apt-get install nodejs
RUN npm install -g pyright@${vPYRIGHT_VERSION}