#!/bin/bash

docker build . -t s21_matrix_test:1.0 -f ./tests/scripts/Dockerfile
docker run s21_matrix_test:1.0
