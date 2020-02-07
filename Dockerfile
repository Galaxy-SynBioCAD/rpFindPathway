FROM brsynth/rpreader-standalone

RUN pip install numpy

COPY rpTool_findpath.py /home/
COPY rpToolServe_findpath.py /home/

COPY results_topass.tar.xz /home/
COPY validation.tar.xz /home/
