(parmmosahn_env) r0715377@set-d-di-09006:~/Documents/parmmosahn$ parmmosahn polish -M test_run/04-curated_models/
Input:  test_run/04-curated_models (4 models)
Output: 05-polished_models
Engine: docker
2026-04-03 14:32:09,967 - parmmosahn.curation - INFO - Polishing 4 models with ModelPolisher (docker)...
2026-04-03 14:32:09,967 - parmmosahn.curation - INFO - Starting local ModelPolisher databases (this may take a moment)...
2026-04-03 14:32:10,049 - parmmosahn.curation - INFO - Waiting 240 seconds for databases to initialize...
2026-04-03 14:36:10,049 - parmmosahn.curation - INFO - Downgrading SBML headers to L3V1 for beta compatibility...
2026-04-03 14:36:10,122 - parmmosahn.curation - INFO - [1/4] Polishing GCA_000227485.1...
2026-04-03 14:36:11,732 - parmmosahn.curation - ERROR - ModelPolisher failed on GCA_000227485.1.
2026-04-03 14:36:11,733 - parmmosahn.curation - ERROR - Docker STDERR:
Cannot connect to the Docker daemon at unix:///home/r0715377/.docker/desktop/docker.sock. Is the docker daemon running?
2026-04-03 14:36:11,733 - parmmosahn.curation - INFO - [2/4] Polishing GCF_000146565.1...
2026-04-03 14:36:11,798 - parmmosahn.curation - ERROR - ModelPolisher failed on GCF_000146565.1.
2026-04-03 14:36:11,798 - parmmosahn.curation - ERROR - Docker STDERR:
permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json?all=1&filters=%7B%22label%22%3A%7B%22com.docker.compose.config-hash%22%3Atrue%2C%22com.docker.compose.project%3Ddata%22%3Atrue%7D%7D": dial unix /var/run/docker.sock: connect: permission denied
2026-04-03 14:36:11,798 - parmmosahn.curation - INFO - [3/4] Polishing GCF_000009045.1...
2026-04-03 14:36:11,856 - parmmosahn.curation - ERROR - ModelPolisher failed on GCF_000009045.1.
2026-04-03 14:36:11,856 - parmmosahn.curation - ERROR - Docker STDERR:
permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json?all=1&filters=%7B%22label%22%3A%7B%22com.docker.compose.config-hash%22%3Atrue%2C%22com.docker.compose.project%3Ddata%22%3Atrue%7D%7D": dial unix /var/run/docker.sock: connect: permission denied
2026-04-03 14:36:11,856 - parmmosahn.curation - INFO - [4/4] Polishing GCF_000209795.2...
2026-04-03 14:36:11,918 - parmmosahn.curation - ERROR - ModelPolisher failed on GCF_000209795.2.
2026-04-03 14:36:11,918 - parmmosahn.curation - ERROR - Docker STDERR:
permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json?all=1&filters=%7B%22label%22%3A%7B%22com.docker.compose.config-hash%22%3Atrue%2C%22com.docker.compose.project%3Ddata%22%3Atrue%7D%7D": dial unix /var/run/docker.sock: connect: permission denied
2026-04-03 14:36:11,918 - parmmosahn.curation - INFO - Shutting down local ModelPolisher databases...
2026-04-03 14:36:11,975 - parmmosahn.curation - WARNING - Polishing finished with 4 failures and 0 successes.

(base) r0715377@set-d-di-09006:~$ docker info
Client: Docker Engine - Community
 Version:    25.0.0
 Context:    default
 Debug Mode: false
 Plugins:
  buildx: Docker Buildx (Docker Inc.)
    Version:  v0.12.0-desktop.2
    Path:     /usr/lib/docker/cli-plugins/docker-buildx
  compose: Docker Compose (Docker Inc.)
    Version:  v2.23.3-desktop.2
    Path:     /usr/lib/docker/cli-plugins/docker-compose
  dev: Docker Dev Environments (Docker Inc.)
    Version:  v0.1.0
    Path:     /usr/lib/docker/cli-plugins/docker-dev
  extension: Manages Docker extensions (Docker Inc.)
    Version:  v0.2.21
    Path:     /usr/lib/docker/cli-plugins/docker-extension
  feedback: Provide feedback, right in your terminal! (Docker Inc.)
    Version:  0.1
    Path:     /usr/lib/docker/cli-plugins/docker-feedback
  init: Creates Docker-related starter files for your project (Docker Inc.)
    Version:  v0.1.0-beta.10
    Path:     /usr/lib/docker/cli-plugins/docker-init
  sbom: View the packaged-based Software Bill Of Materials (SBOM) for an image (Anchore Inc.)
    Version:  0.6.0
    Path:     /usr/lib/docker/cli-plugins/docker-sbom
  scan: Docker Scan (Docker Inc.)
    Version:  v0.26.0
    Path:     /usr/lib/docker/cli-plugins/docker-scan
  scout: Docker Scout (Docker Inc.)
    Version:  v1.2.0
    Path:     /usr/lib/docker/cli-plugins/docker-scout

Server:
ERROR: permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/info": dial unix /var/run/docker.sock: connect: permission denied
errors pretty printing info



