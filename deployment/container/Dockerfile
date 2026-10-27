FROM node:24-bookworm-slim@sha256:a5f4366e834695d5f363b920c104953a535d3a9864d1bdf74566a0ebe01daee5 AS node
FROM docker:cli@sha256:1a4c7cb63513f349bdad01fcc6e0f3f2f67d37b9da86f14dc0d4a0942eecda00 AS dockercli
FROM ghcr.io/astral-sh/uv@sha256:4a1f3666816dafadd20ceb71eac9a7d56d7ef30279f7d64048519ab901ada336 AS uv
FROM python:3.14-slim-bookworm@sha256:48b13b003dda20b16f9442b8475aa05fe21bf6579a8c881db92ffb4d8fd20f83
ENV DEBIAN_FRONTEND=noninteractive PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 UV_NO_PROGRESS=1 NEXT_TELEMETRY_DISABLED=1 MAGNUS_DELIVER=TRUE
COPY --from=node /usr/local/ /usr/local/
COPY --from=dockercli /usr/local/bin/docker /usr/local/bin/docker
COPY --from=uv /uv /usr/local/bin/uv
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates curl tini && rm -rf /var/lib/apt/lists/*
COPY sources/Magnus /opt/magnus
COPY sources/ColliderAgent /opt/collideragent
WORKDIR /opt/magnus/back_end
# Keep the upstream source and lock intact. The control plane needs CPU Torch only.
RUN uv export --frozen --no-dev --no-hashes --no-emit-project --no-emit-workspace --prune torch --prune torchvision > /tmp/backend-requirements.txt && uv pip install --system -r /tmp/backend-requirements.txt
RUN uv pip install --system 'torch==2.9.1+cpu' 'torchvision==0.24.1+cpu' --index-url https://download.pytorch.org/whl/cpu
RUN uv pip install --system --no-deps /opt/magnus/sdks/python
COPY container/build_config.json /opt/magnus/configs/magnus_config.yaml
ENV MAGNUS_CONFIG_PATH=/opt/magnus/configs/magnus_config.yaml
WORKDIR /opt/magnus/front_end
RUN npm ci --no-audit --no-fund && npm run build
COPY container/service.py container/start.sh container/register.py /opt/llp/
COPY container/docker-guard /opt/llp/bin/docker
RUN chmod 755 /opt/llp/bin/docker
ENV PATH=/opt/llp/bin:$PATH
RUN chmod 755 /opt/llp/start.sh && git config --global url.file:///opt/collideragent.insteadOf https://github.com/HET-AGI/ColliderAgent.git && git config --global url.file:///opt/magnus.insteadOf https://github.com/Rise-AGI/magnus.git && git config --global protocol.file.allow always && pip freeze > /opt/llp/python-environment.txt
WORKDIR /opt/magnus/back_end
HEALTHCHECK --interval=15s --timeout=5s --start-period=60s CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8017/health',timeout=3).read()"
ENTRYPOINT ["/usr/bin/tini", "--", "/opt/llp/start.sh"]
