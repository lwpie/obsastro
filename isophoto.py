import astropy.wcs
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.utils.data import download_file
from photutils.isophote import Ellipse, EllipseGeometry

# baseurl = 'https://www.spacetelescope.org/static/projects/fits_liberator/datasets/ngc1068/'
# url1 = baseurl + '791wmos.zip'

# path = download_file(url1)
# hdu = fits.open(path)

# data, header = hdu[0].data, hdu[0].header

hdu = fits.open('fits/10-image-r.fits.fz')[0]
data, header = hdu.data, hdu.header
wcs = astropy.wcs.WCS(header)

# data = ma.masked_equal(data, np.zeros(shape=data.shape))

# g = EllipseGeometry(590, 970, 100, 0.3, 120 / 180 * np.pi)
g = EllipseGeometry(*wcs.wcs.crpix, 100, 0.3, 90 / 180 * np.pi)
g.find_center(data)
ellipse = Ellipse(data, geometry=g)

isolist = ellipse.fit_image(
    integrmode='median', sclip=3.0, nclip=3, fflag=0.3)

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(data, vmin=0, vmax=1200)
ax.set_title("791 wide filter")
ax.set_xlim([300, 900])
ax.set_ylim([700, 1200])

isos = []
for sma in [50., 100., 150., 200., 235.]:
    iso = isolist.get_closest(sma)
    isos.append(iso)
    x, y, = iso.sampled_coordinates()
    plt.plot(x, y, color='w')

plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(isolist.sma**0.25, -2.5*np.log10(isolist.intens))
plt.title("791 wide filter profile")
plt.xlabel('sma**1/4')
plt.ylabel('Magnitude')
plt.gca().invert_yaxis()
plt.show()
