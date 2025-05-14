FROM continuumio/miniconda3:4.6.14
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
ADD HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt /motifs/

RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl libgomp1 build-essential git zlib1g-dev libbz2-dev libcurl4-openssl-dev liblzma-dev libbam-dev libtbb-dev libhdf5-dev libgoogle-perftools-dev libboost-dev libboost-timer-dev libboost-program-options-dev libboost-filesystem-dev libboost-system-dev \
&& conda install -y -c bioconda pandas numpy biopython twobitreader \
&& pip install --upgrade pip \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/*

RUN python -c "from twobitreader import download; download.save_genome('hg38')" \
&& mv hg38.2bit /motifs

RUN git clone https://github.com/thouis/motifscan.git \
&& cd motifscan/liquidator \
&& make CPPFLAGS="-I. -I/usr/local/include -I/usr/include/hdf5/serial"  install \
&& cd ../.. \
&& rm -rf motifscan

ADD scan.py /software/
ADD HOCOMOCOv11_core.meme /motifs/

ENV PATH=/software/:${PATH}
