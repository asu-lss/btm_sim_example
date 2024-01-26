import os
import sys
from glob import glob

import numpy as np
import scipy.linalg as la
import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import random
import h5py
import healpy as hp
from caput import mpiutil, pipeline, config, mpiarray
from draco.core import containers
from drift.core import manager
from drift.core.beamtransfer import matrix_image, matrix_nullspace


def gaussian(x, mu, sigma):
    return (2 * np.pi * sigma ** 2) ** -0.5 * np.exp(-((x - mu) ** 2) / 2 / sigma ** 2)


def fisher_errorbars_from_file(fisher_file, pstype, normalize=True):

    if pstype not in ["unwindowed", "minimum_variance", "uncorrelated"]:
        raise NotImplementedError("pstype %s not implemented!" % pstype)

    # Import stuff from file
    with h5py.File(fisher_file, "r") as f:

        # Fisher matrix
        fisher = f["fisher"][:]
        # Inverse of Fisher matrix (computed by driftscan as
        # la.pinv(self.fisher, rcond=1e-8) )
        fisherinv = f["covariance"][:]

        # Edges of kpar and kperp bins, and numbers of bins
        kpar_edges = f["kpar_bands"][:]
        kperp_edges = f["kperp_bands"][:]
        n_kpar = kpar_edges.shape[0] - 1
        n_kperp = kperp_edges.shape[0] - 1

        # Fisher errorbars for power spectrum bandpowers,
        # for unwindowed power spectrum estimator,
        # packed as [kpar,kperp]
        # (these are just sqrt(diag(fisherinv))
        unwindowed_errs = f["errors"][:].reshape(n_kperp, n_kpar).T

        # kpar and kperp values at bin centers, as 1d arrays or
        # packed as [kpar,kperp]
        kpar_center = f["kpar_center"][:]
        kperp_center = f["kperp_center"][:]
        kpar_center_mesh = kpar_center.reshape(n_kperp, n_kpar).T
        kperp_center_mesh = kperp_center.reshape(n_kperp, n_kpar).T

        # Get bandpowers from file - need to normalize errors by them
        # if unit_bands=False
        band_powers = f["band_power"][:].reshape(n_kperp, n_kpar).T

    if pstype == "unwindowed":

        errs = unwindowed_errs

    elif pstype == "minimum_variance":

        Mab = fisher.sum(axis=1)
        Mab = np.diag(1 / Mab)

        cov = np.dot(Mab, np.dot(fisher, Mab))

        # Fisher errorbars, packed as [kpar,kperp]
        errs = np.sqrt(np.diag(cov)).reshape(n_kperp, n_kpar).T

    elif pstype == "uncorrelated":

        fisher_half = la.cholesky(fisher)
        Mab = la.inv(fisher_half) / fisher_half.sum(axis=1)[:, np.newaxis]

        cov = np.dot(Mab, np.dot(fisher, Mab.T))

        # Fisher errorbars, packed as [kpar,kperp]
        errs = np.sqrt(np.diag(cov)).reshape(n_kperp, n_kpar).T

    # Divide errorbars by fiducial bandpowers. If unit_bands=True,
    # this will just divide by unity, but if unit_bands=False, we'll
    # need to do this division explicitly
    if normalize:
        errs /= band_powers

    return errs, kpar_center_mesh, kperp_center_mesh


def plot_2d_fisher(
    fisher_file, label, pstype, cmin=0.01, cmax=1, cmap="Blues", log=True
):

    if pstype not in ["unwindowed", "minimum_variance", "uncorrelated"]:
        raise NotImplementedError("pstype %s not implemented!" % pstype)

    errs, kpar_mesh, kperp_mesh = fisher_errorbars_from_file(
        fisher_file, pstype, normalize=True
    )

    kpar = kpar_mesh[:, 0]
    kperp = kperp_mesh[0]
    dkpar = kpar[1] - kpar[0]
    dkperp = kperp[1] - kperp[0]

    if log:
        ax = plt.imshow(
            errs,
            origin="lower",
            cmap=cmap,
            vmin=cmin,
            vmax=cmax,
            extent=[
                kperp[0] - 0.5 * dkperp,
                kperp[-1] + 0.5 * dkperp,
                kpar[0] - 0.5 * dkpar,
                kpar[-1] + 0.5 * dkpar,
            ],
            norm=LogNorm(),
        )
    else:
        ax = plt.imshow(
            errs,
            origin="lower",
            cmap=cmap,
            vmin=cmin,
            vmax=cmax,
            extent=[
                kperp[0] - 0.5 * dkperp,
                kperp[-1] + 0.5 * dkperp,
                kpar[0] - 0.5 * dkpar,
                kpar[-1] + 0.5 * dkpar,
            ],
        )
    cbar = plt.colorbar(ax, label=r"$\sigma(P)/P$")
    plt.xticks(
        [
            kperp[0] - 0.5 * dkperp,
            0.5 * (kperp[-1] - kperp[0] + dkperp),
            kperp[-1] + 0.5 * dkperp,
        ]
    )

    # Set plot labels
    #         plt.grid()
    plt.xlabel(r"$k_\perp\; [h\,{\rm Mpc}^{-1}]$")
    plt.ylabel(r"$k_\parallel\; [h\,{\rm Mpc}^{-1}]$")
    plt.title(label)


def plot_fishernormed_measured_ps(
    ps_file,
    fisher_file,
    pstype,
    label,
    n_kpar=0,
    n_kperp=0,
    cmin=-1,
    cmax=1,
    cmap="Blues",
    plot_abs=False,
):

    if pstype not in ["unwindowed", "minimum_variance", "uncorrelated"]:
        raise NotImplementedError("pstype %s not implemented!" % pstype)

    # errs packed as [kpar,kperp]
    errs, kpar_mesh, kperp_mesh = fisher_errorbars_from_file(
        fisher_file, pstype, normalize=False
    )

    kpar = kpar_mesh[:, 0]
    kperp = kperp_mesh[0]
    dkpar = kpar[1] - kpar[0]
    dkperp = kperp[1] - kperp[0]

    # Fetch kpar, kperp, and PS bandpowers from file.
    # Bandpowers are packed as [kperp,kpar], so we transpose
    # to get [kpar,kperp] (axis 1 will be x axis in plot)
    with h5py.File(ps_file, "r") as f:
        #         kpar = f['index_map']['kpar'][:]
        #         kperp = f['index_map']['kperp'][:]
        ps = f["powerspectrum"][:].T

    #     print(errs)
    normed_residuals = ps / errs
    if plot_abs:
        normed_residuals = np.abs(normed_residuals)

    ax = plt.imshow(
        normed_residuals,
        origin="lower",
        cmap=cmap,
        vmin=cmin,
        vmax=cmax,
        extent=[
            kperp[0] - 0.5 * dkperp,
            kperp[-1] + 0.5 * dkperp,
            kpar[0] - 0.5 * dkpar,
            kpar[-1] + 0.5 * dkpar,
        ],
        norm=LogNorm(),
    )
    if plot_abs:
        cbar = plt.colorbar(ax, label=r"$|P - P_{\rm fid}| / \sigma_P$")
    else:
        cbar = plt.colorbar(ax, label=r"$(P - P_{\rm fid}) / \sigma_P$")

    plt.xticks(
        [
            kperp[0] - 0.5 * dkperp,
            0.5 * (kperp[-1] - kperp[0] + dkperp),
            kperp[-1] + 0.5 * dkperp,
        ]
    )

    # Set plot labels
    #         plt.grid()
    plt.xlabel(r"$k_\perp\; [h\,{\rm Mpc}^{-1}]$")
    plt.ylabel(r"$k_\parallel\; [h\,{\rm Mpc}^{-1}]$")
    plt.title(label)
