#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import multiprocessing as mp
import urllib.parse
import os
import sys

import requests

base = 'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/208/2082p145/'
index = 'legacysurvey_dr10_south_coadd_208_2082p145.sha256sum'
base_dir = '2082p145'


def split(url):
    paths = url.split('/')
    base = '/'.join(paths[:-1])
    index = paths[-1]
    base_dir = paths[-2]
    return base, index, base_dir


def join(*iters):
    return '/'.join(iter.rstrip('/') for iter in iters)


def download(url, checksum, base_dir='.'):
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    filename = url.split('/')[-1]
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


def sync(base, index, base_dir='.'):
    recv = requests.get(join(base, index))
    recv.raise_for_status()
    recv.encoding = recv.apparent_encoding
    files = [line.split() for line in recv.text.splitlines()]

    with mp.Pool(8) as pool:
        pool.starmap(download, ((join(base, filename), checksum, base_dir)
                     for checksum, filename in files))


if __name__ == '__main__':
    if len(sys.argv) < 4:
        base, index, base_dir = split(sys.argv[1])
        if len(sys.argv) == 3:
            base_dir = os.path.join(sys.argv[2], base_dir)
        sync(base, index, base_dir)
    elif len(sys.argv) == 4:
        sync(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        sync(base, index, base_dir)
