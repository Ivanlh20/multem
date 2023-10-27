# %%
import h5py
import numpy as np
import os

from PIL import Image
import urllib.request

URL = 'https://upload.wikimedia.org/wikipedia/commons/thumb/8/8b/Denali_Mt_McKinley.jpg/1024px-Denali_Mt_McKinley.jpg'

with urllib.request.urlopen(URL) as url:
    with open('temp.jpg', 'wb') as f:
        f.write(url.read())

img = Image.open('temp.jpg')

# img.show()
# write a image to a hdf5 file in the field "image"
with h5py.File('input.h5', 'w') as f:
    f.create_dataset('image', data=np.array(img)[:,:,0].astype(np.float32))
# %%

# read hdf5 from output.hdf5 in the field "image"
import matplotlib.pyplot as plt
import numpy as np
import h5py
with h5py.File('input.h5', 'r') as f:
    img_in = f['image'][:]
    f.close()

with h5py.File('output.h5', 'r') as f:
    img_out = f['image'][:]
    f.close()
plt.subplot(1,3,1)
plt.imshow(np.array(img_in))
plt.colorbar()
plt.tight_layout()
plt.subplot(1,3,2)
plt.imshow(img_out)
plt.colorbar()
plt.tight_layout()
plt.subplot(1,3,3)
diff= ((img_in-img_out))
plt.imshow(diff)
plt.colorbar()
# %%
