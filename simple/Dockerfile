FROM quay.io/fenicsproject/stable:current

RUN mkdir /home/fenics/pkgs
COPY ./requirements-dev.txt /home/fenics/pkgs/requirements-dev.txt
RUN pip install -r /home/fenics/pkgs/requirements-dev.txt