# Running the EzClermont Webapp

## Docker

### Create the docker image:

```
cd /path/to/ezclermont/webapp/
docker build -t nickp60/ezclermont .
```

### Get an existing docker image:

```
docker pull nickp60/ezclermont
```
### Use the docker image:
```
docker run -p 5000:5000 nickp60/ezclermont
```

### Conda

```

conda env create -f webapp/environment.yml
conda activate clermontwebapp
pip install ezclermont
python clermontwebapp.py

```
