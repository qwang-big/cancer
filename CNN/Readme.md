## resize image
```
for i in {1..20};do d=$((100+$i*10));convert ../img.jpg -resize ${d}% $d.jpg;done
for i in {1..5};do d=$((100-$i*10));convert ../img.jpg -resize ${d}% $d.jpg;done
```
## bash
```
sshfs -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3 192.168.61.8:/prj/ST_PRECISION/OCG/ mnt
```
## get image background
```
from PIL import Image
import numpy as np
import pandas as pd
import sys

f = sys.argv[1]
im = Image.open(f)
a = np.asarray(im)
j = 4 if f[-3:]=='png' else 3
v=np.array([255]*j)-pd.DataFrame(a.reshape(a.shape[0]*a.shape[1],j)).median().values
print(f)
print(",".join(v.astype(np.uint8).astype(np.str)[:3]))
```
## clip image background
```
from PIL import Image
import numpy as np
import pandas as pd
import sys

f = sys.argv[1]
v = sys.argv[2]
v = [float(x) for x in v.split(',')]

im = Image.open(f)
a = np.asarray(im)
d=Image.fromarray(np.clip(a+v,a_min=0,a_max=255).astype(np.uint8))
d.save(f+".jpg")
```
## background percentage
```
from PIL import Image
import numpy as np
import pandas as pd
import sys

f = sys.argv[1]
im = Image.open(f)
gray = im.convert('L')
bw = gray.point(lambda x: 0 if x<220 else 1, 'F')
arr = np.array(np.asarray(bw))
avgBkg = np.average(bw)
bw = gray.point(lambda x: 0 if x<220 else 1, '1')
print('%s %.2f' % (f, avgBkg))
```
