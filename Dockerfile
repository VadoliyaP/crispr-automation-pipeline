FROM python:3.10-slim
RUN useradd -m -u 1000 user
WORKDIR /home/user/app
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential cmake wget ocl-icd-libopencl1 opencl-c-headers pocl-opencl-icd \
    && rm -rf /var/lib/apt/lists/*
RUN wget -q https://github.com/snugel/cas-offinder/archive/refs/heads/master.tar.gz -O cas-offinder.tar.gz && \
    mkdir -p cas-offinder_source && \
    tar -zxf cas-offinder.tar.gz -C cas-offinder_source --strip-components=1 && \
    mkdir -p cas-offinder_source/build && \
    cd cas-offinder_source/build && cmake .. && make && \
    cp cas-offinder /usr/local/bin/cas-offinder && \
    chmod +x /usr/local/bin/cas-offinder && \
    cd /home/user/app && rm -rf cas-offinder.tar.gz cas-offinder_source
COPY --chown=user requirements.txt .
RUN pip install --no-cache-dir --upgrade -r requirements.txt
COPY --chown=user . .
USER user
ENV HOME=/home/user PATH=/home/user/.local/bin:$PATH
EXPOSE 7860
CMD ["python", "app.py"]
