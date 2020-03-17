FROM brsynth/rpreader-standalone:dev

RUN pip install numpy

COPY rpTool_findpath.py /home/
COPY rpToolServe_findpath.py /home/
