#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import astropy.coordinates
import astropy.io.fits
import astropy.stats
import astropy.table
import astropy.units
import astropy.utils.exceptions
import astropy.wcs
import numpy as np
import photutils.aperture
import photutils.background
import photutils.detection
import photutils.segmentation
import scipy.ndimage
from matplotlib import pyplot as plt

import utils

ra, dec = 208.1111458, 14.4908694
band = 'r'

base_dir = 'fits'
figure_dir = 'analyze'


def background(image):
    clipped = astropy.stats.SigmaClip(
        sigma_lower=3, sigma_upper=2, maxiters=10)
    background = photutils.background.Background2D(image, (400, 400), sigma_clip=clipped, bkg_estimator=photutils.background.SExtractorBackground(
        sigma_clip=clipped), bkgrms_estimator=photutils.background.BiweightScaleBackgroundRMS()).background

    mask = np.zeros_like(image, dtype=bool)
    try:
        sources = photutils.detection.DAOStarFinder(
            fwhm=2.5, threshold=np.std(image - background), exclude_border=True).find_stars(image - background)
        mask[sources['ycentroid'].astype(
            int), sources['xcentroid'].astype(int)] = True
        mask = scipy.ndimage.binary_dilation(mask, iterations=5)
    except:
        pass
    background = photutils.background.Background2D(image, (400, 400), mask=mask, sigma_clip=clipped, bkg_estimator=photutils.background.SExtractorBackground(
        sigma_clip=clipped), bkgrms_estimator=photutils.background.BiweightScaleBackgroundRMS()).background

    return background


def magnitude(data):
    return -2.5 * np.log10(data) + 22.5


def aperture(image, invvar, table, wcs):
    points = table[(table['type'] == 'PSF') & (
        table['flux_r'] > 0) & (table['apflux_r'][:, 3] > 0)]

    psf = magnitude(points['flux_r'])
    ap3 = magnitude(points['apflux_r'][:, 3])

    apers = photutils.aperture.SkyCircularAperture(
        astropy.coordinates.SkyCoord(points['ra'], points['dec'], unit='deg'), 1.5 * astropy.units.arcsec)

    phot = photutils.aperture.aperture_photometry(
        image, apers, wcs=wcs, error=np.sqrt(1.0 / invvar), method='subpixel')

    mag = magnitude(phot['aperture_sum'])

    return mag, psf, ap3


def deblend(image, invvar):
    segments = photutils.segmentation.detect_sources(
        image, 5 * np.sqrt(1 / invvar), npixels=10)

    deblend = photutils.segmentation.deblend_sources(
        image, segments, npixels=20, nlevels=8, contrast=0.05)

    return segments, deblend


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        brick = sys.argv[1]
    else:
        brick = ''

    if len(sys.argv) >= 3:
        band = sys.argv[2]

    hdu = astropy.io.fits.open(
        os.path.join(base_dir, f'{brick}image-{band}.fits.fz'))[0]
    image, header = hdu.data, hdu.header
    invvar = astropy.io.fits.open(os.path.join(
        base_dir, f'{brick}invvar-{band}.fits.fz'))[0].data
    table = astropy.table.Table.read(os.path.join(
        base_dir, f'{brick}tractor.fits'), format='fits')
    wcs = astropy.wcs.WCS(header)

    background = background(image)
    utils.plot((background, wcs), filename=os.path.join(
        base_dir, figure_dir, f'{brick}background-{band}.png'))

    mag, psf, ap3 = aperture(image - background, invvar, table, wcs)
    plt.scatter(mag, psf, marker='o', facecolor='none',
                edgecolor='orangered', s=10, label='PSF')
    plt.scatter(mag, ap3, marker='+', color='skyblue', s=20, label='Aper')
    plt.plot([14.5, 29.5], [14.5, 29.5], linestyle='--', color='k')
    plt.legend(loc='best')
    plt.xlabel('Aper')
    plt.ylabel('Tractor')
    utils.finalize(os.path.join(base_dir, figure_dir,
                   f'{brick}aperture-{band}.png'))

    segments, deblend = deblend(image, invvar)
    utils.axis(wcs)
    plt.imshow(deblend, origin='lower',
               interpolation='none', cmap=segments.cmap)
    utils.finalize(os.path.join(base_dir, figure_dir,
                                f'{brick}deblend-{band}.png'))
