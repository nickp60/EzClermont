FROM continuumio/miniconda:4.7.12
MAINTAINER Nick Waters "nickp60@gmail.com"
RUN conda update -y conda && conda install python=3.7
COPY . /
WORKDIR /
RUN pip install -r webapp/requirements.txt
RUN pip install . --use-feature=in-tree-build
ENTRYPOINT ["gunicorn", "-b", "0.0.0.0:5000", "clermontwebapp:app"]
