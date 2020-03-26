FROM brsynth/rpbase:dev

RUN pip install numpy

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY tool_rpFindPathway.py /home/
