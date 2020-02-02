FROM quay.io/fenicsproject/stable:current

RUN mkdir /home/fenics/pkgs
COPY ./requirements-dev.txt /home/fenics/pkgs/requirements-dev.txt
RUN pip install -r /home/fenics/pkgs/requirements-dev.txt

# fenicstools

# (install missing fenicstools deps: https://github.com/mikaem/fenicstools/issues/26)
RUN pip3 install cppimport
RUN apt-get install -y python3-h5py
RUN pip3 install pyvtk

RUN pip3 install git+https://github.com/mikaem/fenicstools