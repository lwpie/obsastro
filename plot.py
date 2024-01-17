#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import astropy.io.fits
import astropy.stats
import astropy.visualization
import numpy as np
import photutils.background
from matplotlib import pyplot as plt

import utils

bands = ['g', 'r', 'i', 'z']
brick = 'image'
base_dir = 'fits'
figure_dir = 'figures'
mapping = True


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


def plot(hdu, filename=None):
    utils.plot(hdu, filename=filename)


def hist(hdu, filename=None):
    image, header = hdu.data, hdu.header
    plt.hist(image.flatten(), bins=100, range=(-0.05, 0.05))
    utils.finalize(filename)


def color(hdus, bound=None, mapping=True, bands=None, filename=None):
    if (not mapping) and (not bands):
        raise ValueError('bands must be specified when mapping is False')
    if not bands:
        bands = list(hdus.keys())
    datas = [hdus[band].data for band in bands]

    if mapping:
        image = dr_rgb(datas, bands)
    else:
        image = astropy.visualization.make_lupton_rgb(
            *datas, stretch=0.5, Q=10)

    utils.plot(astropy.io.fits.PrimaryHDU(image, next(
        iter(hdus.values())).header), finish=not bound, filename=filename)

    if bound:
        data = np.mean(datas, axis=0)
        upper = np.nanpercentile(data.flatten(), bound * 100)
        data[data > upper] = np.nan
        scale = np.arcsinh(data)

        vmin, vmax = astropy.visualization.ZScaleInterval(
            contrast=bound).get_limits(data)
        plt.imshow(scale, vmin=vmin, vmax=vmax, cmap='gray', interpolation='none',
                   origin='lower', clim=[np.min(data), upper], alpha=0.9)

        utils.finalize(filename)


def background(hdu, filename=None):
    image, header = hdu.data, hdu.header
    wcs = astropy.wcs.WCS(header)
    utils.axis(wcs)

    clipped = astropy.stats.SigmaClip(
        sigma_lower=3.0, sigma_upper=2.0, maxiters=10)
    background = photutils.background.Background2D(image, (400, 400), sigma_clip=clipped, bkg_estimator=photutils.background.SExtractorBackground(
        sigma_clip=clipped), bkgrms_estimator=photutils.background.BiweightScaleBackgroundRMS())
    image = background.background

    plt.imshow(image, origin='lower', interpolation='none', vmin=np.nanpercentile(
        image.flatten(), 1), vmax=np.nanpercentile(image.flatten(), 99))
    background.plot_meshes(outlines=True, marker='.', color='cyan', alpha=0.3)

    utils.finalize(filename)


if __name__ == '__main__':
    # plt.style.use('seaborn-v0_8')

    if len(sys.argv) > 1:
        brick = sys.argv[1]
        if len(sys.argv) > 2:
            bands = sys.argv[2].split(',')
            mapping = False

    hdus = {band: astropy.io.fits.open(os.path.join(
        base_dir, f'{brick}-{band}.fits.fz'))[0] for band in bands}

    path = os.path.join(base_dir, figure_dir)

    hist(next(iter(hdus.values())), os.path.join(
        path, f'{brick}-hist.png'))
    plot(next(iter(hdus.values())), os.path.join(
        path, f'{brick}-plot.png'))

    color(hdus, mapping=mapping, bands=bands,
          filename=os.path.join(path, f'{brick}-color.png'))
    color(hdus, 0.8, mapping, bands, os.path.join(
        path, f'{brick}-bound.png'))

    background(next(iter(hdus.values())), os.path.join(
        path, f'{brick}-background.png'))
