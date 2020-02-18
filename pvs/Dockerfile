# Container for offscreen software rending
#
# see EGL based container for cases using the
# nvidia runtime

ARG BASE_IMAGE=kitware/paraviewweb:pvw-v5.6.0-osmesa
FROM ${BASE_IMAGE}

RUN mkdir /opt/paraview/share/paraview-5.6/data

ENTRYPOINT [ "/opt/paraview/bin/pvpython", "/opt/paraview/share/paraview-5.6/web/visualizer/server/pvw-visualizer.py" ]

CMD [ "--port", "8080", "--content", "/opt/paraview/share/paraview-5.6/web/visualizer/www/",\
      "--data", "/opt/paraview/share/paraview-5.6/data" ]