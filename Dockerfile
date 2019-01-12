FROM python:3.6-jessie

# install packages
RUN apt-get update && apt-get install -y \
    vim

# install pytest for testing
RUN pip install pytest

# set up working directory
COPY . /phyllite
WORKDIR /phyllite

CMD /bin/bash
