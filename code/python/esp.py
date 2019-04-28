# Maximilian N. Diefenbach <maximilian.diefenbach@tum.de>
# Body Magnetic Resonance Research Group, http://www.bmrrgroup.de
# October 2018

import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from numpy.linalg import norm
from scipy.interpolate import splrep, splev
from scipy.ndimage.morphology import binary_erosion
from utils import *
from skimage.measure import compare_ssim, compare_mse, compare_nrmse
from mpl_toolkits.axes_grid1 import make_axes_locatable


def index_coords(arr, origin=None, dtype=float):
    """Creates coordinate arrays I, J, K for the indicies in a numpy array "arr".
    "origin" defaults to the center of the image2. Specify origin=(0,0,0)
    to set the origin to the lower left corner of the image2."""
    ni, nj, nk = arr.shape
    I, J, K = np.mgrid[:ni, :nj, :nk]
    if origin is None:
        origin = np.array(arr.shape) // 2
    I -= origin[0]
    J -= origin[1]
    K -= origin[2]
    return I.astype(dtype), J.astype(dtype), K.astype(dtype)


def shell_mask(arr, radius, origin=None):
    I, J, K = index_coords(arr, origin)
    R = np.sqrt(I**2 + J**2 + K**2)
    mask = ((radius-1) <= R) & (R < radius)
    return mask


def cone_mask(arr, angle_rad, direction, origin=None):
    I, J, K = index_coords(arr, origin)
    R = np.sqrt(I**2 + J**2 + K**2)
    B = direction[0] * I + direction[1] * J + direction[2] * K
    return np.abs(B) == np.round(R * np.cos(angle_rad))


def cone_distance(arr, direction, origin=None):
    I, J, K = index_coords(arr, origin)
    R = np.sqrt(I**2 + J**2 + K**2)
    B = direction[0] * I + direction[1] * J + direction[2] * K
    Theta = np.arccos(np.divide(B, R, out=np.zeros_like(B), where=R!=0))
    magic_angle_rad = np.arccos(1/np.sqrt(3))
    return np.array(
        [np.abs(np.round(R * np.sin(magic_angle_rad + Theta))),
         np.abs(np.round(R * np.sin(magic_angle_rad - Theta)))]).min(0)


def esp_shell(img1, img2):
    img1F = ifftshift(fftn(img1))
    err = img1 - img2
    errF = ifftshift(fftn(err))

    I, J, K = index_coords(img1, origin=None)
    R = np.sqrt(I**2 + J**2 + K**2)

    radii = np.arange(1, np.max(img2.shape) // 2 + 1)
    fr = np.zeros_like(radii)
    yr = np.zeros_like(radii)
    for i, radius in enumerate(radii):
        mask = ((radius-1) <= R) & (R < radius) # shell_mask(img1, radius)
        fr[i] = norm(errF[mask])
        yr[i] = norm(img1F[mask])

    return radii, esp(fr, yr)


def esp_cone(img1, img2, direction):
    img1F = ifftshift(fftn(img1))
    err = img1 - img2
    errF = ifftshift(fftn(err))

    I, J, K = index_coords(img1, origin=None)
    R = np.sqrt(I**2 + J**2 + K**2)
    B = direction[0] * I + direction[1] * J + direction[2] * K
    Theta = np.arccos(np.divide(B, R, out=np.zeros_like(B), where=R!=0))
    magic_angle_rad = np.arccos(1/np.sqrt(3))

    angles = np.arange(0, 91)
    fr = np.zeros_like(angles)
    yr = np.zeros_like(angles)
    deg2rad = np.pi / 180
    for i, angle in enumerate(angles):
        angle_rad = angle * deg2rad
        mask = (np.abs(B) == np.round(R * np.sin(np.pi/2 - angle_rad)))
        fr[i] = norm(errF[mask])
        yr[i] = norm(img1F[mask])

    return angles, esp(fr, yr)


def esp_zerocone_distance(img1, img2, direction):
    img1F = ifftshift(fftn(img1))
    err = img1 - img2
    errF = ifftshift(fftn(err))

    zc_dist = cone_distance(img1, direction)
    dist = np.arange(zc_dist.min(), zc_dist.max(), 1)
    fr = np.zeros_like(dist)
    yr = np.zeros_like(dist)
    for i, d in enumerate(dist):
        mask = zc_dist == d
        fr[i] = norm(errF[mask])
        yr[i] = norm(img1F[mask])

    return dist, esp(fr, yr)


def esp(fr, yr, smoothing=1):
    return np.sqrt(splines(fr, smoothing) /
                   splines(yr, smoothing))


def splines(arr1d, smoothing):
    x = np.arange(len(arr1d))
    tck = splrep(x, arr1d, s=smoothing)
    return splev(x, tck, der=0)


if __name__ == '__main__':
    # load data
    f_img1 = 'cosmos_flipped.nii.gz'
    f_img2 = 'chi33_flipped.nii.gz'

    image1 = sitk.ReadImage(f_img1)
    image2 = sitk.ReadImage(f_img2)

    img1 = sitk.GetArrayFromImage(image1).transpose((1, 2, 0))
    img2 = sitk.GetArrayFromImage(image2).transpose((1, 2, 0))

    close_all()
    plt.close('all')

    # compute error spectrum plots (esp)
    radii, Er_shell = esp_shell(img1, img2)
    plt.figure()
    plt.plot(radii[1:], Er_shell[1:], 'o-')
    plt.xlabel('radius [px]')
    plt.ylabel('nrmse')

    B0dir = [0, 0, 1]
    angles_deg, Er_cone = esp_cone(img1, img2, B0dir)
    plt.figure()
    plt.plot(angles_deg, Er_cone, 'o-')
    plt.xlabel('angle(B0, theta) [deg]')
    plt.ylabel('nrmse')

    dist, Er_zc_dist = esp_zerocone_distance(img1, img2, B0dir)
    plt.figure()
    plt.plot(dist, Er_zc_dist, 'o-')
    plt.xlabel('zero-cone distance [px]')
    plt.ylabel('nrmse')

    radius_maxerror = radii[Er_shell.argmax()]
    radius2_maxerror = radii[1:][Er_shell[1:].argmax()]
    angles_deg_maxerror = angles_deg[Er_cone.argmax()]
    dist_maxerror = Er_zc_dist.argmax()

    print(f'r_max = {radius_maxerror}')
    print(f'r2_max = {radius2_maxerror}')
    print(f'alpha_max = {angles_deg_maxerror}')
    print(f'dist_max = {dist_maxerror}')

    ssim = compare_ssim(img1, img2)
    mse = compare_mse(img1, img2)
    nrmse = compare_nrmse(img1, img2)
    print(f'ssim = {ssim},\nmse = {mse},\nnrmse = {nrmse}')

    zc_dist = cone_distance(img1, B0dir)
    show_arr3d(zc_dist == dist_maxerror)

    # figures

    # compute FFT
    img1F = ifftshift(fftn(img1))
    img2F = ifftshift(fftn(img2))

    fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(2, 3, figsize=(17,20))

    im0 = ax0.imshow(img1[:, :, img1.shape[2]//2])
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im0, cax=cax, orientation='vertical')
    ax0.set_title('COSMOS')

    im1 = ax1.imshow(img2[:, :, img2.shape[2]//2])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    ax1.set_title('STI 33')

    im2 = ax2.imshow((img1 - img2)[:, :, img2.shape[2]//2])
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical')
    ax2.set_title('COSMOS - STI 33')

    im3 = ax3.imshow(np.log(np.abs(img1F))[:, :, img1F.shape[2]//2])
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3, cax=cax, orientation='vertical')
    ax3.set_title('log |FFT COSMOS|')

    im4 = ax4.imshow(np.log(np.abs(img2F[:, :, img2F.shape[2]//2])))
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im4, cax=cax, orientation='vertical')
    ax4.set_title('log |FFT STI 33|')

    im5 = ax5.imshow(np.log(np.abs(fftn(img1-img2)))[:, :, img2F.shape[2]//2])
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im5, cax=cax, orientation='vertical')
    ax5.set_title('log | |FFT COSMOS| - |FFT STI 33| |')

    for ax in [ax0, ax1, ax2, ax3, ax4, ax5]:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.tight_layout()

    show_arr3d(np.log(np.abs(fftn(img1-img2))))
