#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import astropy.coordinates
import astropy.io.fits
import astropy.nddata
import astropy.wcs
import numpy as np
import reproject
import reproject.mosaicking
from matplotlib import pyplot as plt

channel = 'image'
bands = ['r', 'g', 'i', 'z']
bricks = ['2079p145', '2082p145']
base_dir = 'fits'
ra, dec = 208.1111458, 14.4908694


def axis(wcs):
    plt.figure(figsize=(10, 5))
    plt.subplot(projection=wcs)
    plt.xlabel('R.A [deg]')
    plt.ylabel('Dec [deg]')
    # plt.grid(color='white', ls='solid')


def finalize(filename=None):
    path = os.path.dirname(filename)
    if not os.path.exists(path):
        os.makedirs(path)

    if filename:
        plt.savefig(filename, dpi=300)
    else:
        plt.show()

    plt.close()


def merge(filenames):
    data = [astropy.io.fits.open(filename)[1] for filename in filenames]
    if len(data) == 1:
        return data[0]
    wcs, shape = reproject.mosaicking.find_optimal_celestial_wcs(data)
    array, footprint = reproject.mosaicking.reproject_and_coadd(
        data, wcs, shape_out=shape, reproject_function=reproject.reproject_interp)
    return astropy.io.fits.PrimaryHDU(array, header=wcs.to_header())


def crop(hdu, center, scale):
    image, header = hdu.data, hdu.header
    wcs = astropy.wcs.WCS(header)

    wcs.wcs.crval = center
    wcs.wcs.pc = [[1 / scale, 0], [0, 1 / scale]]
    image, footprint = reproject.reproject_interp(
        hdu, wcs)

    return astropy.io.fits.PrimaryHDU(image, header=wcs.to_header())

    # center = astropy.coordinates.SkyCoord(*center, unit='deg')
    # size = np.array(image.shape) / scale
    # hdu = astropy.nddata.Cutout2D(
    #     image, center, size, wcs=wcs)

    # return astropy.io.fits.PrimaryHDU(hdu.data, header=hdu.wcs.to_header())


def plot(hdu, finish=True, filename=None):
    image, header = hdu.data, hdu.header
    wcs = astropy.wcs.WCS(header, naxis=2)

    axis(wcs)
    plt.imshow(image, vmin=np.nanpercentile(image.flatten(), 10), vmax=np.nanpercentile(
        image.flatten(), 99),  origin='lower', interpolation='none')

    if finish:
        finalize(filename)


if __name__ == '__main__':
    # plt.style.use('seaborn-v0_8')

    if len(sys.argv) >= 2:
        channel = sys.argv[1]
        if len(sys.argv) >= 3:
            bricks = sys.argv[3].split(',')
    for band in bands:
        filenames = [os.path.join(
            base_dir, brick, f'{channel}-{band}.fits.fz') for brick in bricks]
        hdu = crop(merge(filenames), (ra, dec), 10)
        plot(hdu, filename=os.path.join(
            base_dir, f'{channel}-{band}.png'))
        hdu.writeto(os.path.join(
            base_dir, f'{channel}-{band}.fits.fz'), overwrite=True)
