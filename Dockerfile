FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for the IKMB VEP pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/vep-2.0/bin:$PATH

RUN apt-get -y update && apt-get -y install make wget unzip build-essential libyaml-dev libssl-dev zlib1g-dev
RUN cd /opt && wget https://cache.ruby-lang.org/pub/ruby/3.2/ruby-3.2.2.tar.gz && tar -xvf ruby-3.2.2.tar.gz && cd ruby-3.2.2 && ./configure && make -j2 install
RUN gem install fast_excel
