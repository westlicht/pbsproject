#!/usr/bin/env python

import json
import random

random.seed(1337041)

scale = 10

particleRadius = scale * 0.05

w = scale * 20
h = scale * 5
d = scale * 10
border = scale * 1
water = scale * 1
waterheight = scale * 3
waterborder = scale * 0.01

bcount = 100
bminsize = scale * 0.1
bmaxsize = scale * 0.5
bminheight = scale * 0.5
bmaxheight = scale * 3

boxes = []

# buidings
for i in range(bcount):
    bw = random.uniform(bminsize, bmaxsize)
    bh = random.uniform(bminheight, bmaxheight)
    bd = random.uniform(bminsize, bmaxsize)
    x = random.uniform(-w/2+bw/2+border,w/2-bw/2-border-water*3)
    z = random.uniform(-d/2+bd/2+border,d/2-bd/2-border)
    boxes.append({
        "bounds" : [[x-bw,0,z-bd],[x+bw,bh,z+bd]],
        "type" : "boundary"
    })

# water
boxes.append({
    "bounds" : [[w/2-waterborder-water,0,-d/2+waterborder],[w/2-waterborder,waterheight,d/2-waterborder]],
    "type" : "fluid"
})

settings = {
    "particleRadius" : particleRadius,
    "viscosity" : 0.0000001 / (scale*scale*scale),
}

scene = {
    "camera" : {
        "fov" : 50,
        "position" : [0,w/1.5,d],
        "target" : [0,0,0],
        "up" : [0,0,-1],
        "near" : 0.1,
        "far" : 2*w
    },
    "world" : {
        "bounds" : [[-w/2,0,-d/2],[w/2,h,d/2]]
    },
    "boxes" : boxes
}

data = { "settings" : settings, "scene" : scene }

print json.dumps(data)

with open('tsunami.json', 'w') as outfile:
    json.dump(data, outfile)
