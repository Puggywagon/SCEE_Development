#!/usr/bin/python3

import yaml


# this is a test code to read a YAML file



yamlfile = 'example1.yml'
with open(yamlfile, 'r') as f:
    tmp = yaml.safe_load_all(f)
    for t in tmp:
        print(t)
