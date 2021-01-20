
############################ specify base image ############################
FROM ubuntu:18.04


############################ copy stuff from the host machine into the container image ############################
COPY ./gemse_analysis_docker_image /home/gemse_analysis_infrastructure


############################ specify shell ############################
SHELL ["/bin/bash", "-c"]


############################ execute commands from within the container to install ROOT, BAT and Moritz' scripts ############################
RUN apt-get update \
    && apt-get --assume-yes upgrade \
    ############################ installing jupyter notebook ############################
    && apt-get --assume-yes install python \
    && apt-get --assume-yes install python3 \
    && apt-get --assume-yes install python3-pip \
    && pip3 install numpy \
    && pip3 install jupyter
    #&& apt-get --assume-yes install jupyter-notebook 
    ############################ installing ROOT ############################
    #&& apt-get --assume-yes install make \
    #&& apt-get --assume-yes install cmake \
    #&& apt-get --assume-yes install gcc \
    #&& apt-get --assume-yes install g++ \
    #&& apt-get --assume-yes install binutils \
    #&& apt-get --assume-yes install libx11-dev \
    #&& apt-get --assume-yes install libxpm-dev \
    #&& apt-get --assume-yes install libxft-dev \
    #&& apt-get --assume-yes install libxext-dev \
    #&& apt-get --assume-yes install libssl-dev \
    #&& apt-get --assume-yes install vim \
    #&& /bin/bash -c "source /home/gemse_analysis_infrastructure/root/root_v6.22.02.Linux-ubuntu18-x86_64-gcc7.5/root/bin/thisroot.sh" \
    #&& echo "export PATH=$PATH:~/.local/bin" >>~/.bashrc \
    #&& echo "source /home/gemse_analysis_infrastructure/root/root_v6.22.02.Linux-ubuntu18-x86_64-gcc7.5/root/bin/thisroot.sh" >>~/.bashrc \
    ############################# installing BAT ############################
    #&& cd /home/gemse_analysis_infrastructure/bat/BAT-0.9.4.1 \
    #&& source /home/gemse_analysis_infrastructure/root/root_v6.22.02.Linux-ubuntu18-x86_64-gcc7.5/root/bin/thisroot.sh \
    #&& ./configure --prefix=/home/gemse_analysis_infrastructure/bat_install \
    #&& make \
    #&& make install \
    #&& echo 'BATPREFIX="/home/gemse_analysis_infrastructure/bat_install"' >>~/.bashrc \
    #&& echo 'export PATH="$BATPREFIX/bin:$PATH"' >>~/.bashrc \
    #&& echo 'export LD_LIBRARY_PATH="$BATPREFIX/lib:$LD_LIBRARY_PATH"' >>~/.bashrc \
    #&& echo 'export CPATH="$BATPREFIX/include:$CPATH"' >>~/.bashrc \
    #&& echo 'export PKG_CONFIG_PATH="$BATPREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"' >>~/.bashrc \
    #&& source ~/.bashrc \
    #&& BATPREFIX="/home/gemse_analysis_infrastructure/bat_install" \
    #&& export PATH="$BATPREFIX/bin:$PATH" \
    #&& export LD_LIBRARY_PATH="$BATPREFIX/lib:$LD_LIBRARY_PATH" \
    #&& export CPATH="$BATPREFIX/include:$CPATH" \
    #&& export PKG_CONFIG_PATH="$BATPREFIX/lib/pkgconfig:$PKG_CONFIG_PATH" \
    #&& echo $SHELL \
    #&& echo $PATH \
    ############################# installing Moritz' scripts ############################
    #&& cd /home/gemse_analysis_infrastructure/gemse_root_scripts \
    #&& make \
    #&& cd /home/gemse_analysis_infrastructure/gemse_analysis \
    #&& make
    ############################# installing Python and Jupyter Notebook ############################



