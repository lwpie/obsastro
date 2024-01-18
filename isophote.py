#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import astropy.coordinates
import astropy.io.fits
import astropy.wcs
import numpy as np
import photutils.isophote
import scipy.ndimage
from matplotlib import pyplot as plt

import utils

masx = 208.1111458, 14.4908694
ic = 208.1068042, 14.4886816
band = 'r'

base_dir = 'fits'
figure_dir = 'isophote'


def mask(image, center, sma, angle, **kwargs):
    g = photutils.isophote.EllipseGeometry(
        *center, sma, 0., angle * np.pi / 180)
    g.find_center(image)
    ellipse = photutils.isophote.Ellipse(image, geometry=g)
    isolist = ellipse.fit_image(**kwargs)

    model = photutils.isophote.build_ellipse_model(image.shape, isolist)
    return model


def isophote(image, center, sma, angle, **kwargs):
    g = photutils.isophote.EllipseGeometry(
        *center, sma, 0., angle * np.pi / 180)
    g.find_center(image)

    ellipse = photutils.isophote.Ellipse(image, geometry=g)
    return ellipse.fit_image(**kwargs)


if __name__ == '__main__':
    # plt.style.use('seaborn-v0_8')

    if len(sys.argv) >= 2:
        brick = sys.argv[1]
    else:
        brick = ''

    if len(sys.argv) >= 3:
        band = sys.argv[2]

    hdu = astropy.io.fits.open(
        os.path.join(base_dir, f'{brick}image-{band}.fits.fz'))[0]
    image, header = np.ma.masked_equal(
        hdu.data, np.zeros(shape=hdu.data.shape)), hdu.header
    wcs = astropy.wcs.WCS(header)

    model = mask(image, utils.to_pix(*ic, wcs), 100, 30, maxsma=300)
    image -= model
    utils.plot((image - model, wcs), filename=os.path.join(
        base_dir, figure_dir, f'{brick}masked-{band}.png'))

    isolist = isophote(image, utils.to_pix(*masx, wcs),
                       400, 100, maxsma=500, step=0.3)
    utils.plot((image, wcs), finish=False)
    for sma in np.arange(50, isolist.sma.max(), 50):
        iso = isolist.get_closest(sma)
        x, y = iso.sampled_coordinates()
        plt.plot(x, y, color='w')
    utils.finalize(os.path.join(
        base_dir, figure_dir, f'{brick}isophote-{band}.png'))

    area = np.pi * isolist.sma ** 2 * (1 - isolist.eps)
    flux = np.append(area[0], np.diff(area)) * isolist.intens

    plt.figure(figsize=(8, 4))
    plt.scatter(isolist.sma[1:] ** 0.25,
                utils.magnitude(isolist.intens)[1:], label='Point')
    plt.scatter(isolist.sma[1:] ** 0.25,
                utils.magnitude(np.cumsum(flux))[1:], label='Cumulative')
    plt.xlabel('sma ** 1/4')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.legend(loc='best')
    utils.finalize(os.path.join(
        base_dir, figure_dir, f'{brick}profile-{band}.png'))

    lines = list()
    for angle in np.arange(0, 360, 30):
        rotated = scipy.ndimage.rotate(image, angle, reshape=False)
        y, x = utils.to_pix(*masx, wcs)
        line = rotated[x, y: y + 1000]
        lines.append(line := line / line.max())
        # plt.plot(np.diff(line), label=f'{angle} deg')

    plt.figure(figsize=(10, 5))
    for line in lines:
        plt.plot(line, label=f'{angle} deg')
    plt.legend(loc='upper right')
    utils.finalize(os.path.join(
        base_dir, figure_dir, f'{brick}grad-{band}.png'))

    plt.figure(figsize=(10, 5))
    for line in lines:
        plt.plot(np.diff(line), label=f'{angle} deg')
    plt.legend(loc='lower right')
    utils.finalize(os.path.join(
        base_dir, figure_dir, f'{brick}gradiff-{band}.png'))

    # plt.figure(figsize=(10, 5))
    # plt.figure(1)

    # plt.subplot(221)
    # plt.errorbar(isolist.sma, isolist.eps,
    #              yerr=isolist.ellip_err, fmt='o', markersize=4)
    # plt.xlabel('Semimajor axis length')
    # plt.ylabel('Ellipticity')

    # plt.subplot(222)
    # plt.errorbar(isolist.sma, isolist.pa/np.pi*180.,
    #              yerr=isolist.pa_err/np.pi * 80., fmt='o', markersize=4)
    # plt.xlabel('Semimajor axis length')
    # plt.ylabel('PA (deg)')

    # plt.subplot(223)
    # plt.errorbar(isolist.sma, isolist.x0,
    #              yerr=isolist.x0_err, fmt='o', markersize=4)
    # plt.xlabel('Semimajor axis length')
    # plt.ylabel('X0')

    # plt.subplot(224)
    # plt.errorbar(isolist.sma, isolist.y0,
    #              yerr=isolist.y0_err, fmt='o', markersize=4)
    # plt.xlabel('Semimajor axis length')
    # plt.ylabel('Y0')

    # plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10,
    #                     right=0.95, hspace=0.35, wspace=0.35)

    # utils.finalize(os.path.join(
    #     base_dir, figure_dir, f'{brick}ellipse-{band}.png'))
