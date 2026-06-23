FROM python:3.10-slim

WORKDIR /code

# Install basic system requirements if needed, but no heavy drivers
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Setup user permissions for Hugging Face Spaces
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH

WORKDIR $HOME/app

COPY --chown=user . $HOME/app

RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "app.py"]
