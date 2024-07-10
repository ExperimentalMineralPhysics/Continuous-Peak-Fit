#!/usr/bin/env python

import sys
from PIL import Image
from PIL.ExifTags import TAGS
from cpf.XRD_FitPattern import logger

Im = sys.argv[1]

logger.info(" ".join(map(str, [("Image File: ", Im)])))

with Image.open(Im) as img:
    meta_dict = {TAGS[key]: img.tag[key] for key in img.tag.iterkeys()}

for (tag, value) in Image.open(Im)._getexif().iteritems():
    logger.info(" ".join(map(str, [("%s = %s" % (TAGS.get(tag), value))])))
