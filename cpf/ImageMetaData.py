#!/usr/bin/env python
import os
import sys
from PIL import Image
from PIL.ExifTags import TAGS

Im = sys.argv[1]

print 'Image File: ', Im

with Image.open(Im) as img:
    meta_dict = {TAGS[key] : img.tag[key] for key in img.tag.iterkeys()}

for (tag,value) in Image.open(Im)._getexif().iteritems():
        print '%s = %s' % (TAGS.get(tag), value)