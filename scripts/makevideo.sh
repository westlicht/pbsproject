#!/bin/sh

ffmpeg -i frame%04d.png -vcodec libx264 -crf 5 fluid.mp4

