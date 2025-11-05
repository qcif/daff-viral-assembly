#!/usr/bin/env bash

TAG=neoformit/daff-ont-assembly:latest

docker build -t $TAG . && docker push $TAG
