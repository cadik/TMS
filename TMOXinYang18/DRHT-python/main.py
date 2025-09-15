#import img_io as iio
#import hdr2ldr as h2l
#import ldr2hdr as l2h

import sys
import numpy as np

print("RUN", file=sys.stderr)
print(sys.path)
if len(sys.argv) != 4:
  sys.exit("WrongNumberOfArrrrrgs", file=sys.stderr)
imgPath = str(sys.argv[1])
imgW = int(sys.argv[2])
imgH = int(sys.argv[3])
print(imgPath, file=sys.stderr)
print(imgW, file=sys.stderr)
print(imgH, file=sys.stderr)
#loadedimg = iio.readLDR(imgPath,[imgW,imgH])
#iio.writeLDR(loadedimg,imgPath+".jpg")