#!/bin/bash

name='pi@192.168.0.113:/home/pi/'
rsync -vr --exclude "fundamental*" --exclude "CMakeFiles" --exclude "CMakeCache*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/python_shared/* ${name}/research/pi/
