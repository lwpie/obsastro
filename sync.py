#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import functools
import hashlib
import multiprocessing as mp
import os
import sys

import requests

# import urllib.parse


base = r'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/{channel}/{brick:.3}/{brick}/'
index = r'legacysurvey_dr10_south_{channel}_{brick:.3}_{brick}.sha256sum'
brick = '2082p145'
base_dir = 'fits'


def split(url):
    paths = url.split('/')
    base = '/'.join(paths[:-1])
    index = paths[-1]
    base_dir = paths[-2]
    return base, index, base_dir


def join(*iters):
    return '/'.join(iter.rstrip('/') for iter in iters)


def download(url, checksum, base_dir='.', fix=''):
    if not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    filename = os.path.basename(url).replace(fix, '')
    path = os.path.join(base_dir, filename)
    print(path, file=sys.stderr)

    if os.path.exists(path):
        with open(path, 'rb') as f:
            data = f.read()
    else:
        recv = requests.get(url)
        recv.raise_for_status()
        with open(path, 'wb') as f:
            f.write(data := recv.content)

    if hashlib.sha256(data).hexdigest() != checksum:
        raise ValueError(f'Checksum mismatch: {filename}')


def sync(base, index, base_dir='.', sync=True):
    recv = requests.get(join(base, index))
    recv.raise_for_status()
    recv.encoding = recv.apparent_encoding

    files = dict()
    for line in recv.text.splitlines():
        checksum, filename = line.split()
        files[join(base, filename.strip('*'))] = checksum
    prefix = os.path.commonprefix(list(map(os.path.basename, files)))

    if sync:
        with mp.Pool(8) as pool:
            pool.starmap(download, ((*file, base_dir, prefix)
                         for file in files.items()))

    return files


if __name__ == '__main__':
    if len(sys.argv) == 2:
        brick = sys.argv[1]
        if len(sys.argv) == 3:
            base_dir = sys.argv[2]

    base = functools.partial(base.format, brick=brick)
    index = functools.partial(index.format, brick=brick)
    path = os.path.join(base_dir, brick)

    sync(base(channel='coadd'), index(channel='coadd'), path)

    base = base(channel='tractor').replace(brick, '')
    index = index(channel='tractor').replace(f'_{brick}', '')
    filename = join(base, f'tractor-{brick}.fits')
    download(filename, sync(base, index, sync=False)
             [filename], path, f'-{brick}')
