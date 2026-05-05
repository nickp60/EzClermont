FROM python:3.12-slim

WORKDIR /app

COPY webapp/requirements.txt .
RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

COPY . .

RUN pip install --no-cache-dir .

RUN useradd -m appuser
USER appuser

EXPOSE 5000

ENTRYPOINT ["gunicorn", "-b", "0.0.0.0:5000", "webapp.clermontwebapp:app"]
