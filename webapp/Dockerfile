FROM ubuntu:latest
MAINTAINER Nick Waters "nickp60@gmail.com"
RUN apt-get update -y
RUN apt-get install curl -y

RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3/
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda3/bin:${PATH}
RUN conda update -y conda
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt
RUN pip install ezclermont
ENTRYPOINT ["python"]
CMD ["clermontwebapp.py"]