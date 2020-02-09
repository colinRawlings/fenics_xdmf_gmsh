FROM ceciledc/fenics_mixed_dimensional:27-01-20 

ARG vGMSH_VERSION=4.5.0
ARG vPYRIGHT_VERSION=1.1.19

RUN apt-get update

RUN apt-get install -y libgfortran3 libglu1-mesa libxcursor1 libxinerama1
RUN mkdir /home/fenics/apps/
RUN curl -o /home/fenics/apps/gmsh.tgz http://gmsh.info/bin/Linux/gmsh-${vGMSH_VERSION}-Linux64-sdk.tgz
RUN cd /home/fenics/apps && tar -xvzf /home/fenics/apps/gmsh.tgz
ENV PYTHONPATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/lib:${PYTHONPATH}"
ENV PATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/bin:${PATH}"

RUN pip install meshio
RUN pip install lxml
RUN pip install --no-binary=h5py h5py

RUN mkdir /home/fenics/pkgs
COPY ./requirements-dev.txt /home/fenics/requirements-dev.txt
RUN pip install -r /home/fenics/requirements-dev.txt

RUN curl -sL https://deb.nodesource.com/setup_12.x | sudo -E bash -
RUN apt-get install nodejs
RUN npm install -g pyright@${vPYRIGHT_VERSION}

# fenicstools

# (install missing fenicstools deps: https://github.com/mikaem/fenicstools/issues/26)
RUN pip3 install cppimport
RUN apt-get install -y python3-h5py
RUN pip3 install pyvtk

RUN pip3 install git+https://github.com/mikaem/fenicstools

# (trigger import to get fenicstools to do its compilation)
RUN python3 -c "import fenicstools"


CMD [ "/bin/bash" ]