FROM brsynth/rpbase

RUN pip install numpy pandas sklearn

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpFindPathway.py /home/
