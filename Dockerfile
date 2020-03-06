FROM ubuntu:latest
MAINTAINER stayaside
RUN apt-get update && apt-get install -y apt-transport-https python3.6
RUN mkdir /data
ADD lec2_hw_sample1.fsa /data/lec2_hw_sample1.fsa
ADD lec2_hw_sample3.fsa /data/lec2_hw_sample3.fsa
ADD exercise1.py /exercise1.py
ENV TARGET_FILE /data/lec2_hw_sample1.fsa
ENV TARGET_FILE_2 /data/lec2_hw_sample3.fsa
CMD cat $TARGET_FILE | python3.6 exercise1.py && cat $TARGET_FILE_2 | python3.6 exercise1.py

