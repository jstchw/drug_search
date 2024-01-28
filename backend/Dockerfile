FROM python:3.10-slim

WORKDIR /usr/src/app

COPY requirements.txt /usr/src/app/

RUN apt-get update && apt-get install -y \
    build-essential \
    libxrender1 \
    libxext6 \
    libfontconfig1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir -r requirements.txt

COPY . /usr/src/app

EXPOSE 16000

CMD ["python3", "-m", "flask", "run", "--host=0.0.0.0", "-p 16000"]