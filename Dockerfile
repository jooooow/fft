FROM alpine:3.20.3

RUN apk add --no-cache make

RUN apk add g++

COPY fft_1d.cpp /fft_1d.cpp
COPY Makefile ./Makefile

CMD make all
