#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
import os
import sys

import utils
import astropy.visualization
import astropy.io.fits
import numpy as np

bands = ['g', 'r', 'i', 'z']
brick = 'merged'
base_dir = 'fits'
mapping = True


def denoise(image):
    return image
    vmin, vmax = np.nanpercentile(image.flatten(), (1, 99))
    image = np.clip(image, vmin, vmax)
    return image


def plot(hdu, filename=None):
    return utils.plot(hdu, filename=filename)


def color(hdus, mapping=True, bands=None, filename=None):
    if (not mapping) and (not bands):
        raise ValueError('bands must be specified when mapping is False')
    if not bands:
        bands = list(hdus.keys())
    datas = [denoise(hdus[band].data) for band in bands]
    if mapping:
        image = dr_rgb(datas, bands)
    else:
        image = astropy.visualization.make_lupton_rgb(
            *datas, stretch=0.5, Q=10)
    # flatten = image.sum(axis=2)
    # vmin = np.percentile(flatten, 10)
    # index = flatten < vmin
    # image[index] = (0, 0, 0)
    return utils.plot(astropy.io.fits.PrimaryHDU(image, next(iter(hdus.values())).header), filename=filename)


def dr_rgb(imgs, bands):
    m = 0.03
    Q = 20
    f = 1.5
    scales = dict(g=6.0, r=3.4, i=3.0, z=2.2)

    I = 0
    for img, band in zip(imgs, bands):
        I += np.maximum(0, img * scales[band] * f + m)
    I /= len(bands)
    fI = np.arcsinh(Q * I) / np.sqrt(Q)
    I += (I == 0.) * 1e-6
    I = fI / I
    H, W = I.shape
    rgb = np.zeros((H, W, 3), np.float32)

    vec = dict(
        g=(0.,   0.,  0.75),
        r=(0.,   0.5, 0.25),
        i=(0.25, 0.5, 0.),
        z=(0.75, 0.,  0.))

    for img, band in zip(imgs, bands):
        scale = scales[band] * f
        rf, gf, bf = vec[band]
        v = np.clip((img * scale + m) * I, 0, 1)
        if rf != 0.:
            rgb[:, :, 0] += rf * v
        if gf != 0.:
            rgb[:, :, 1] += gf * v
        if bf != 0.:
            rgb[:, :, 2] += bf * v
    return rgb


if __name__ == '__main__':
    # plt.style.use('seaborn-v0_8')

    if len(sys.argv) > 1:
        brick = sys.argv[1]
        if len(sys.argv) > 2:
            filename = sys.argv[2]
        if len(sys.argv) > 3:
            bands = sys.argv[3].split(',')
            mapping = False

    hdus = {band: astropy.io.fits.open(os.path.join(
        base_dir, f'{brick}-{band}.fits.fz'))[0] for band in bands}
    color(hdus, mapping, bands, filename)
