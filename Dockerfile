FROM python:3.13

RUN mkdir /code

WORKDIR /code

COPY requirements.txt /code/

RUN python -m pip install -U pip setuptools

RUN python -m pip install -r requirements.txt

RUN playwright install-deps
RUN playwright install

RUN apt purge -y fonts-wqy-zenhei
RUN apt install -y fonts-noto

ENV PYDEVD_DISABLE_FILE_VALIDATION=1
