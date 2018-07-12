# Running the EzClermont Webapp

## Docker

Create the docker image:

```
cd /path/to/ezclermont/webapp/
docker build -t ezclermont .
docker run docker run -p 5000:5000 ezclermont
```


### Virtualenv

run within a python virtual 
```
venv ~/ezenv/ -p python3
source ~/ezenv/bin/activate
pip install flask ezclermont biopython
cd /path/to/ezclermont/webapp/
python clermontwebapp.py

```