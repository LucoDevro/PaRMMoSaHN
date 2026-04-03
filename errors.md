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




(parmmosahn_env) r0715377@set-d-di-09006:~/Documents/parmmosahn$ parmmosahn run -g genomes \
> -m M9_glc_aerobic.csv \
> -l subtilus \
> -o test_run/
Input:  genomes
Output: test_run
Label:  subtilus

Steps 1-2: Building pangenome and pan-model...
2026-04-03 11:16:18,565 - parmmosahn.pipeline - INFO - === STEP 1: Building Pangenome ===
2026-04-03 11:16:18,565 - parmmosahn.pangenome - INFO - Checking for problematic characters in GenBank files in genomes...
2026-04-03 11:16:18,590 - parmmosahn.pangenome - INFO - Generating annotations list at test_run/00-pan_results/01-pangenome/annotations.list...
2026-04-03 11:16:18,590 - parmmosahn.pangenome - INFO - Running PPanGGOLiN annotate...
2026-04-03 11:16:19,587 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:19 utils.py:l239 INFO	Command: /home/r0715377/miniconda3/envs/parmmosahn_env/bin/ppanggolin annotate -f --anno test_run/00-pan_results/01-pangenome/annotations.list -o test_run/00-pan_results/01-pangenome -c 14
2026-04-03 11:16:19,588 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:19 utils.py:l242 INFO	PPanGGOLiN version: 2.2.6
2026-04-03 11:16:19,589 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:19 utils.py:l977 INFO	3 parameters have a non-default value.
2026-04-03 11:16:19,589 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:19 annotate.py:l1533 INFO	Reading annotations.list the list of genome files ...
2026-04-03 11:16:19,928 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:19 annotate.py:l1579 INFO	gene identifiers used in the provided annotation files were unique, PPanGGOLiN will use them.
2026-04-03 11:16:19,929 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:19 writeBinaries.py:l921 INFO	Writing genome annotations...
2026-04-03 11:16:21,321 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:21 writeBinaries.py:l932 INFO	writing the protein coding gene dna sequences in pangenome...
2026-04-03 11:16:22,010 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:22 writeMetadata.py:l411 INFO	Writing genomes metadata in pangenome
2026-04-03 11:16:22,022 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:22 writeMetadata.py:l428 INFO	Writing contigs metadata in pangenome
2026-04-03 11:16:22,037 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:22 writeBinaries.py:l995 INFO	Done writing the pangenome. It is in file : test_run/00-pan_results/01-pangenome/pangenome.h5
2026-04-03 11:16:22,243 - parmmosahn.pangenome - INFO - Running PPanGGOLiN cluster...
2026-04-03 11:16:23,253 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 utils.py:l239 INFO	Command: /home/r0715377/miniconda3/envs/parmmosahn_env/bin/ppanggolin cluster -f -p test_run/00-pan_results/01-pangenome/pangenome.h5 --identity 0.8 --coverage 0.6 -c 14
2026-04-03 11:16:23,254 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 utils.py:l242 INFO	PPanGGOLiN version: 2.2.6
2026-04-03 11:16:23,255 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 utils.py:l977 INFO	3 parameters have a non-default value.
2026-04-03 11:16:23,255 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 readBinaries.py:l123 INFO	Getting the current pangenome status
2026-04-03 11:16:23,257 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 readBinaries.py:l346 INFO	Extracting and writing CDS sequences from a /home/r0715377/Documents/parmmosahn/test_run/00-pan_results/01-pangenome/pangenome.h5 file to a fasta file...
2026-04-03 11:16:23,465 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 cluster.py:l478 INFO	Clustering all of the genes sequences...
2026-04-03 11:16:23,466 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 writeSequences.py:l117 INFO	Creating sequence database...
2026-04-03 11:16:23,490 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:23 cluster.py:l136 INFO	Clustering sequences...
2026-04-03 11:16:24,760 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:24 cluster.py:l163 INFO	Extracting cluster representatives...
2026-04-03 11:16:24,776 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:24 cluster.py:l185 INFO	Writing gene to family information
2026-04-03 11:16:24,794 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:24 cluster.py:l488 INFO	Associating fragments to their original gene family...
2026-04-03 11:16:24,794 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:24 writeSequences.py:l117 INFO	Creating sequence database...
2026-04-03 11:16:24,804 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:24 cluster.py:l249 INFO	Aligning cluster representatives...
2026-04-03 11:16:26,657 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 cluster.py:l274 INFO	Extracting alignments...
2026-04-03 11:16:26,676 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 cluster.py:l332 INFO	Starting with 5397 families
2026-04-03 11:16:26,682 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 cluster.py:l372 INFO	Ending with 5328 gene families
2026-04-03 11:16:26,684 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 cluster.py:l385 INFO	Adding protein sequences to the gene families
2026-04-03 11:16:26,728 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 cluster.py:l402 INFO	Adding 16278 genes to the gene families
2026-04-03 11:16:26,751 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 cluster.py:l861 INFO	Done with the clustering
2026-04-03 11:16:26,751 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 writeBinaries.py:l939 INFO	Writing gene families and gene associations in pangenome...
2026-04-03 11:16:26,782 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 writeBinaries.py:l943 INFO	Writing gene families information in pangenome...
2026-04-03 11:16:26,858 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:26 writeBinaries.py:l743 INFO	Updating annotations with fragment information
2026-04-03 11:16:27,007 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:27 writeBinaries.py:l995 INFO	Done writing the pangenome. It is in file : /home/r0715377/Documents/parmmosahn/test_run/00-pan_results/01-pangenome/pangenome.h5
2026-04-03 11:16:27,179 - parmmosahn.pangenome - INFO - Extracting pangenome statistics...
2026-04-03 11:16:28,179 - parmmosahn.utils - INFO - [ppanggolin] Content:
2026-04-03 11:16:28,179 - parmmosahn.utils - INFO - [ppanggolin] Genes: 16278
2026-04-03 11:16:28,179 - parmmosahn.utils - INFO - [ppanggolin] Genomes: 4
2026-04-03 11:16:28,179 - parmmosahn.utils - INFO - [ppanggolin] Families: 5328
2026-04-03 11:16:28,339 - parmmosahn.pangenome - INFO - Calculated soft_core threshold: 0.5 for 4 genomes.
2026-04-03 11:16:28,339 - parmmosahn.pangenome - INFO - Extracting soft-core protein families...
2026-04-03 11:16:29,342 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:29 utils.py:l239 INFO	Command: /home/r0715377/miniconda3/envs/parmmosahn_env/bin/ppanggolin fasta -f -p test_run/00-pan_results/01-pangenome/pangenome.h5 -o test_run/00-pan_results/01-pangenome --prot_families softcore --soft_core 0.5
2026-04-03 11:16:29,343 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:29 utils.py:l242 INFO	PPanGGOLiN version: 2.2.6
2026-04-03 11:16:29,344 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:29 utils.py:l977 INFO	4 parameters have a non-default value.
2026-04-03 11:16:29,344 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:29 readBinaries.py:l123 INFO	Getting the current pangenome status
2026-04-03 11:16:29,346 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:29 writeSequences.py:l484 INFO	Writing the representative protein sequences of the gene families by reading the pangenome file directly.
2026-04-03 11:16:29,415 - parmmosahn.utils - INFO - [ppanggolin] 2026-04-03 11:16:29 readBinaries.py:l747 INFO	Done writing the representative amino acid sequences of the gene families:'test_run/00-pan_results/01-pangenome/softcore_protein_families.faa'
2026-04-03 11:16:29,582 - parmmosahn.pangenome - INFO - Pangenome construction complete! File saved to: test_run/00-pan_results/01-pangenome/subtilus.faa
2026-04-03 11:16:29,582 - parmmosahn.pipeline - INFO - === STEP 2: Building Pan-model ===
2026-04-03 11:16:29,582 - parmmosahn.models - INFO - Cleaning formatting for subtilus.faa...
2026-04-03 11:16:29,585 - parmmosahn.models - INFO - Running gapseq find. This may take a while depending on the pangenome size...
2026-04-03 13:18:40,969 - parmmosahn.models - INFO - Pan-model generation complete! Results saved in /home/r0715377/Documents/parmmosahn/test_run/00-pan_results/02-panmodel
2026-04-03 13:18:40,970 - parmmosahn.pipeline - INFO - Pan-tasks complete! Pan-model saved to test_run/00-pan_results/02-panmodel
Step 3: Deriving strain-specific models...
2026-04-03 13:18:40,970 - parmmosahn.models - INFO - Preparing input genomes for DIAMOND mapping...
2026-04-03 13:18:41,183 - parmmosahn.models - INFO - Cleaning formatting for GCA_000227485.1.faa...
2026-04-03 13:18:41,466 - parmmosahn.models - INFO - Cleaning formatting for GCF_000209795.2.faa...
2026-04-03 13:18:41,778 - parmmosahn.models - INFO - Cleaning formatting for GCF_000146565.1.faa...
2026-04-03 13:18:42,125 - parmmosahn.models - INFO - Cleaning formatting for GCF_000009045.1.faa...
2026-04-03 13:18:42,129 - parmmosahn.models - INFO - Building Diamond database from pangenome...
2026-04-03 13:18:42,170 - parmmosahn.models - INFO - [1/4] --- Processing strain: GCA_000227485.1 ---
2026-04-03 13:24:34,926 - parmmosahn.models - INFO - [2/4] --- Processing strain: GCF_000209795.2 ---
2026-04-03 13:30:23,764 - parmmosahn.models - INFO - [3/4] --- Processing strain: GCF_000146565.1 ---
2026-04-03 13:36:21,781 - parmmosahn.models - INFO - [4/4] --- Processing strain: GCF_000009045.1 ---
2026-04-03 13:42:27,507 - parmmosahn.models - INFO - All strain models successfully generated in /home/r0715377/Documents/parmmosahn/test_run/03-strain_models/models
Steps 4-7: Curation, polishing, annotation, gather, final evaluation...
2026-04-03 13:42:27,509 - parmmosahn.pipeline - INFO - === STEP 4a: Generating curation template (memote on drafts) ===
2026-04-03 13:42:27,509 - parmmosahn.curation - INFO - Running automated memote evaluation on 4 draft models...
2026-04-03 13:42:27,510 - parmmosahn.curation - INFO -   [1/4] Evaluating GCA_000227485.1...
2026-04-03 13:42:27,510 - parmmosahn.curation - INFO - Running memote snapshot on GCA_000227485.1.xml...
2026-04-03 13:50:07,544 - parmmosahn.curation - INFO - Memote snapshot completed for GCA_000227485.1.xml.
2026-04-03 13:50:07,544 - parmmosahn.curation - INFO -   [2/4] Evaluating GCF_000009045.1...
2026-04-03 13:50:07,544 - parmmosahn.curation - INFO - Running memote snapshot on GCF_000009045.1.xml...
2026-04-03 13:57:49,983 - parmmosahn.curation - INFO - Memote snapshot completed for GCF_000009045.1.xml.
2026-04-03 13:57:49,983 - parmmosahn.curation - INFO -   [3/4] Evaluating GCF_000146565.1...
2026-04-03 13:57:49,983 - parmmosahn.curation - INFO - Running memote snapshot on GCF_000146565.1.xml...
2026-04-03 14:05:40,469 - parmmosahn.curation - INFO - Memote snapshot completed for GCF_000146565.1.xml.
2026-04-03 14:05:40,469 - parmmosahn.curation - INFO -   [4/4] Evaluating GCF_000209795.2...
2026-04-03 14:05:40,469 - parmmosahn.curation - INFO - Running memote snapshot on GCF_000209795.2.xml...
2026-04-03 14:06:40,816 - parmmosahn.curation - INFO - Memote snapshot completed for GCF_000209795.2.xml.
2026-04-03 14:06:40,816 - parmmosahn.curation - INFO - Extracting problem balances and duplicates from memote reports...
2026-04-03 14:06:40,919 - parmmosahn.curation - INFO - Writing curation template to curation_template.xlsx...
2026-04-03 14:06:43,324 - parmmosahn.pipeline - INFO - Curation template generated: test_run/curation_template.xlsx
2026-04-03 14:06:43,324 - parmmosahn.pipeline - INFO - === STEP 4: Applying Automatic Curations ===
2026-04-03 14:06:57,778 - parmmosahn.pipeline - INFO - === STEP 5: ModelPolisher (docker) ===
2026-04-03 14:06:57,778 - parmmosahn.curation - INFO - Polishing 4 models with ModelPolisher (docker)...
2026-04-03 14:06:57,779 - parmmosahn.curation - INFO - Starting local ModelPolisher databases (this may take a moment)...
2026-04-03 14:06:57,851 - parmmosahn.curation - INFO - Waiting 240 seconds for databases to initialize...
2026-04-03 14:10:57,852 - parmmosahn.curation - INFO - Downgrading SBML headers to L3V1 for beta compatibility...
2026-04-03 14:10:57,916 - parmmosahn.curation - INFO - [1/4] Polishing GCA_000227485.1...
2026-04-03 14:11:06,221 - parmmosahn.curation - ERROR - ModelPolisher failed on GCA_000227485.1.
2026-04-03 14:11:06,221 - parmmosahn.curation - ERROR - Docker STDERR:
Cannot connect to the Docker daemon at unix:///home/r0715377/.docker/desktop/docker.sock. Is the docker daemon running?
2026-04-03 14:11:06,221 - parmmosahn.curation - INFO - [2/4] Polishing GCF_000146565.1...
2026-04-03 14:11:06,289 - parmmosahn.curation - ERROR - ModelPolisher failed on GCF_000146565.1.
2026-04-03 14:11:06,290 - parmmosahn.curation - ERROR - Docker STDERR:
permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json?all=1&filters=%7B%22label%22%3A%7B%22com.docker.compose.config-hash%22%3Atrue%2C%22com.docker.compose.project%3Ddata%22%3Atrue%7D%7D": dial unix /var/run/docker.sock: connect: permission denied
2026-04-03 14:11:06,290 - parmmosahn.curation - INFO - [3/4] Polishing GCF_000009045.1...
2026-04-03 14:11:06,352 - parmmosahn.curation - ERROR - ModelPolisher failed on GCF_000009045.1.
2026-04-03 14:11:06,352 - parmmosahn.curation - ERROR - Docker STDERR:
permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json?all=1&filters=%7B%22label%22%3A%7B%22com.docker.compose.config-hash%22%3Atrue%2C%22com.docker.compose.project%3Ddata%22%3Atrue%7D%7D": dial unix /var/run/docker.sock: connect: permission denied
2026-04-03 14:11:06,352 - parmmosahn.curation - INFO - [4/4] Polishing GCF_000209795.2...
2026-04-03 14:11:06,413 - parmmosahn.curation - ERROR - ModelPolisher failed on GCF_000209795.2.
2026-04-03 14:11:06,413 - parmmosahn.curation - ERROR - Docker STDERR:
permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json?all=1&filters=%7B%22label%22%3A%7B%22com.docker.compose.config-hash%22%3Atrue%2C%22com.docker.compose.project%3Ddata%22%3Atrue%7D%7D": dial unix /var/run/docker.sock: connect: permission denied
2026-04-03 14:11:06,413 - parmmosahn.curation - INFO - Shutting down local ModelPolisher databases...
2026-04-03 14:11:06,474 - parmmosahn.curation - WARNING - Polishing finished with 4 failures and 0 successes.
2026-04-03 14:11:06,474 - parmmosahn.pipeline - INFO - === STEP 5b: Annotation Enrichment ===
2026-04-03 14:11:06,474 - parmmosahn.annotation - WARNING - No XML models found in /home/r0715377/Documents/parmmosahn/test_run/05-polished/models — skipping annotation.
2026-04-03 14:11:06,474 - parmmosahn.pipeline - INFO - === STEP 6: Gathering Final Results ===
2026-04-03 14:11:06,474 - parmmosahn.results - WARNING - No .xml models found in /home/r0715377/Documents/parmmosahn/test_run/05b-annotated to gather.
2026-04-03 14:11:06,474 - parmmosahn.curation - WARNING - No final XML models found to evaluate.
2026-04-03 14:11:06,474 - parmmosahn.pipeline - INFO - Pipeline complete! All results saved to test_run/06-final_results

Pipeline complete! (2h 54m 47s total)
Results: test_run



