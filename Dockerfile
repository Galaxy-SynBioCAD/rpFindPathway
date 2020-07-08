FROM brsynth/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \
                        libxext6  \
        libxrender-dev \
         && conda install -y -c rdkit rdkit

RUN pip install numpy pandas sklearn

COPY rpFindPathway.py /home/
COPY rpFindPathwayServe.py /home/
COPY run_rpFindPathway.py /home/
