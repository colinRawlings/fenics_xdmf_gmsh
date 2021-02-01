ARG vFENICS_CONTAINER=ceciledc/fenics_mixed_dimensional:05-02-20

FROM ${vFENICS_CONTAINER}

ARG vGMSH_VERSION=4.5.0
ARG vPYRIGHT_VERSION=1.1.19

# gmsh

RUN apt-get update && apt-get install -y libgfortran3 libglu1-mesa libxcursor1 libxinerama1
RUN mkdir /home/fenics/apps/
RUN curl -o /home/fenics/apps/gmsh.tgz http://gmsh.info/bin/Linux/gmsh-${vGMSH_VERSION}-Linux64-sdk.tgz
RUN cd /home/fenics/apps && tar -xvzf /home/fenics/apps/gmsh.tgz
ENV PYTHONPATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/lib:${PYTHONPATH}"
ENV PATH="/home/fenics/apps/gmsh-${vGMSH_VERSION}-Linux64-sdk/bin:${PATH}"

# pyright

RUN curl -sL https://deb.nodesource.com/setup_12.x | sudo -E bash -
RUN apt-get install nodejs
RUN npm install -g pyright@${vPYRIGHT_VERSION}

# fenicstools

# (install missing fenicstools deps: https://github.com/mikaem/fenicstools/issues/26)
RUN pip3 install cppimport==18.11.8
RUN apt-get install -y python3-h5py
RUN pip3 install pyvtk==0.5.18

RUN pip3 install git+https://github.com/mikaem/fenicstools@2019.1

# (trigger import to get fenicstools to do its compilation)
RUN python3 -c "import fenicstools"

# paraview

RUN apt install -y paraview \ 
                   x11-apps
COPY ./docker/fenics/ParaView-UserSettings.json /root/.config/ParaView/ParaView-UserSettings.json

# jupyterlab + notebook

RUN pip install --upgrade jupyterlab
ENV JUPYTERLAB_DIR="$HOME/.local/share/jupyter/lab"
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager
RUN jupyter labextension install jupyter-matplotlib
RUN jupyter nbextension enable --py widgetsnbextension

RUN pip install jupytext --upgrade

# homemade mesh conversion

COPY ./docker/fenics/requirements.txt /home/fenics/requirements.txt
RUN pip install -r /home/fenics/requirements.txt

COPY ./docker/fenics/requirements-dev.txt /home/fenics/requirements-dev.txt
RUN pip install -r /home/fenics/requirements-dev.txt

COPY ./pkgs /home/fenics/pkgs
RUN pip install /home/fenics/pkgs/fenics_utils
RUN pip install /home/fenics/pkgs/paraview_utils
# for development it is likely preferrable to overwrite these with an editably installed copy

ENV MPLBACKEND="TkAgg"

CMD [ "/bin/bash" ]