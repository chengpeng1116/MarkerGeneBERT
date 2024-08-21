FROM lfoppiano/grobid:0.8.0

# 克隆GitHub仓库
RUN apt-get update --allow-releaseinfo-change
RUN apt-get install -y git --fix-missing

# Install dependencies, including libbz2-dev for bz2 support
RUN apt-get update --fix-missing
RUN apt-get install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev  libxml2-dev libcurl4-openssl-dev 
RUN apt-get install -y wget unzip libbz2-dev  procps vim

# Download and install Python with bz2 support
RUN apt-get install -y r-base  && \
	wget https://mirrors.huaweicloud.com/python/3.9.9/Python-3.9.9.tgz && \
    tar -xvf Python-3.9.9.tgz && \
    cd Python-3.9.9 && \
    ./configure --enable-optimizations --with-bz2 && \
    make altinstall 

RUN update-alternatives --install /usr/bin/python3 python3 /usr/local/bin/python3.9 1

# Install pip
RUN apt-get install -y python3-pip --fix-missing

WORKDIR /opt/spacy_model
RUN 
	pip3 install tqdm pandas dataclasses tk spacy scipdf  scispacy  inflect sci-hub -i https://pypi.tuna.tsinghua.edu.cn/simple && \
	pip3 install git+https://github.com/titipata/scipdf_parser -i https://pypi.tuna.tsinghua.edu.cn/simple
	wget https://github.com/explosion/spacy-models/releases/download/en_core_web_sm-3.7.1/en_core_web_sm-3.7.1.tar.gz && \
	wget https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_scibert-0.5.4.tar.gz && \
	wget https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_jnlpba_md-0.5.4.tar.gz && \
	wget https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_craft_md-0.5.4.tar.gz && \
	wget https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz && \
	wget https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bc5cdr_md-0.5.4.tar.gz && \
	pip3 install /opt/spacy_model/en_core_web_sm-3.7.1.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple && \
	pip3 install /opt/spacy_model/en_core_sci_scibert-0.5.4.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple && \
	pip3 install /opt/spacy_model/en_ner_jnlpba_md-0.5.4.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple && \
	pip3 install /opt/spacy_model/en_ner_craft_md-0.5.4.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple && \
	pip3 install /opt/spacy_model/en_ner_bionlp13cg_md-0.5.4.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple && \
	pip3 install /opt/spacy_model/en_ner_bc5cdr_md-0.5.4.tar.gz -i https://pypi.tuna.tsinghua.edu.cn/simple  && \
	R -e "install.packages('xml2',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" && \
	R -e "install.packages('curl',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" && \
	R -e "install.packages('RISmed',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" && \
    R -e "install.packages('optparse',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"  && \
	R -e "install.packages('easyPubMed',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"  && \
	R -e "install.packages('tidypmc',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"  && \
	R -e "install.packages('tidyverse',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"  && \
	R -e "install.packages('europepmc',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"  && \
    R -e "install.packages('rentrez',dependencies=TRUE, repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" 


WORKDIR /opt//docker_file
RUN git clone  https://github.com/chengpeng1116/MarkerGeneBERT.git

RUN chmod +x /opt/docker_file/MarkerGeneBERT/script/start.sh
# 启动应用
CMD ["/opt//docker_file/MarkerGeneBERT/script//start.sh"]









