Bootstrap: docker
From: continuumio/miniconda3

%labels
    MAINTAINER Carlos Guzman <cag104@eng.ucsd.edu>
    DESCRIPTION Container image containing all requirements for the adapted bennerlab/smrnaseq pipeline
    VERSION 0.1dev

%files
    environment.yml /

%environment
	PATH=/opt/conda/envs/smrnaseq-0.1dev/bin:$PATH
	export PATH

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
