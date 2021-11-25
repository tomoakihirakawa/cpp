#!/bin/bash

name='pi@10.0.1.7:/home/pi/'
rsync -vr --exclude "fundamental*" --exclude "CMakeFiles" --exclude "CMakeCache*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/python_shared/* ${name}/research/pi/
