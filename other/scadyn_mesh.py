#!/usr/bin/env python3
# Copyright 2019 Joonas Herranen (joonas.herranen@iki.fi)
#
# Mesh utilities for scadyn: HDF5 I/O, Gaussian sphere/ellipsoid generation, layered
# permittivity, and visualization. CLI subcommands or GUI (default).

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits import mplot3d
from numpy.random import normal, seed
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform
from scipy.special import factorial, lpmv, spherical_in

# ---------------------------------------------------------------------------
# HDF5 I/O
# ---------------------------------------------------------------------------


def _vertices_from_coord(coord: np.ndarray) -> np.ndarray:
    coord = np.asarray(coord, dtype=np.float64)
    if coord.ndim != 2:
        raise ValueError(f"coord must be 2-D, got shape {coord.shape}")
    if coord.shape[0] == 3 and coord.shape[1] != 3:
        return coord.T
    if coord.shape[1] == 3:
        return coord
    raise ValueError(f"unrecognized coord shape {coord.shape}")


def _tets_from_etopol(etopol: np.ndarray) -> np.ndarray:
    etopol = np.asarray(etopol, dtype=np.int64)
    if etopol.ndim != 2:
        raise ValueError(f"etopol must be 2-D, got shape {etopol.shape}")
    if etopol.shape[0] == 4:
        tets = etopol.T
    elif etopol.shape[1] >= 4:
        tets = etopol[:, :4]
    else:
        raise ValueError(f"unrecognized etopol shape {etopol.shape}")
    return tets - 1


def read_h5(path: str) -> dict:
    if not path.endswith(".h5"):
        path = path + ".h5"
    with h5py.File(path, "r") as f:
        vertices = _vertices_from_coord(f["coord"][:])
        tets = _tets_from_etopol(f["etopol"][:])
        out: Dict[str, Any] = {
            "vertices": vertices,
            "tets": tets,
            "path": path,
        }
        if "param_r" in f:
            out["param_r"] = np.asarray(f["param_r"][:], dtype=np.float64)
        if "param_i" in f:
            out["param_i"] = np.asarray(f["param_i"][:], dtype=np.float64)
    return out


def write_h5(
    path: str,
    vertices: np.ndarray,
    tets: np.ndarray,
    param_r: np.ndarray,
    param_i: np.ndarray,
) -> str:
    if not path.endswith(".h5"):
        path = path + ".h5"
    vertices = np.asarray(vertices, dtype=np.float64)
    tets = np.asarray(tets, dtype=np.int32)
    if tets.ndim != 2 or tets.shape[1] != 4:
        raise ValueError(f"tets must be (M, 4), got {tets.shape}")
    n_tet = tets.shape[0]
    param_r = np.broadcast_to(np.asarray(param_r, dtype=np.float64), (n_tet,))
    param_i = np.broadcast_to(np.asarray(param_i, dtype=np.float64), (n_tet,))
    with h5py.File(path, "w") as f:
        f.create_dataset("coord", data=vertices)
        f.create_dataset("etopol", data=tets + 1)
        f.create_dataset("param_r", data=param_r)
        f.create_dataset("param_i", data=param_i)
    return path


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def boundary_faces(tets: np.ndarray) -> np.ndarray:
    tets = np.asarray(tets, dtype=np.int64)
    faces = np.vstack(
        [
            tets[:, [0, 1, 2]],
            tets[:, [0, 1, 3]],
            tets[:, [0, 2, 3]],
            tets[:, [1, 2, 3]],
        ]
    )
    faces = np.sort(faces, axis=1)
    unique, inverse = np.unique(faces, axis=0, return_inverse=True)
    counts = np.bincount(inverse)
    return unique[counts == 1]


def tetra_volume(vertices: np.ndarray, tets: np.ndarray) -> float:
    vol = 0.0
    for tet in tets:
        p0, p1, p2, p3 = vertices[tet]
        vol += abs(np.dot(p0 - p3, np.cross(p1 - p3, p2 - p3))) / 6.0
    return vol


def effective_radius(vertices: np.ndarray, tets: np.ndarray) -> float:
    return (3.0 * tetra_volume(vertices, tets) / (4.0 * np.pi)) ** (1.0 / 3.0)


def validate_surface_mesh(
    vertices: np.ndarray, faces: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return contiguous surface arrays safe for TetGen."""
    vertices = np.ascontiguousarray(vertices, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError(f"vertices must be (N, 3), got {vertices.shape}")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError(f"faces must be (M, 3), got {faces.shape}")
    if vertices.shape[0] < 4 or faces.shape[0] < 4:
        raise ValueError("surface mesh too small for tetrahedralization")
    if not np.isfinite(vertices).all():
        raise ValueError("surface mesh has non-finite vertex coordinates")
    n = vertices.shape[0]
    if faces.min() < 0 or faces.max() >= n:
        raise ValueError(
            f"face indices out of range (n_vertices={n}, max index={faces.max()})"
        )
    return vertices, faces


def tetrahedralize(
    vertices: np.ndarray,
    faces: np.ndarray,
    refinement: float,
) -> Tuple[np.ndarray, np.ndarray]:
    import pyvista as pv
    import tetgen

    vertices, faces = validate_surface_mesh(vertices, faces)
    if refinement <= 0:
        raise ValueError("refinement (max tet volume) must be positive")
    pv_faces = np.ascontiguousarray(
        np.hstack([np.full((faces.shape[0], 1), 3, dtype=np.int64), faces])
    )
    surface = pv.PolyData(vertices, pv_faces)
    if not surface.is_all_triangles:
        surface = surface.triangulate()
    tgen = tetgen.TetGen(surface)
    nodes, elems, _, _ = tgen.tetrahedralize(
        order=1,
        maxvolume=float(refinement),
        minratio=1.5,
        mindihedral=20.0,
    )
    return np.asarray(nodes, dtype=np.float64), np.asarray(elems, dtype=np.int64)


def face_vectors(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    return vertices[faces]


def format_mesh_info(data: dict) -> str:
    v, t = data["vertices"], data["tets"]
    vol = tetra_volume(v, t)
    a_eff = effective_radius(v, t)
    lines = [
        f"File:      {data['path']}",
        f"Nodes:     {v.shape[0]}",
        f"Tets:      {t.shape[0]}",
        f"Volume:    {vol:.6g}",
        f"a_eff:     {a_eff:.6g}",
    ]
    if "param_r" in data:
        pr = data["param_r"]
        lines.append(
            f"param_r:   min={pr.min():.6g} max={pr.max():.6g} "
            f"unique={np.unique(pr).size}"
        )
    if "param_i" in data:
        pi = data["param_i"]
        lines.append(
            f"param_i:   min={pi.min():.6g} max={pi.max():.6g} "
            f"unique={np.unique(pi).size}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plotting (matplotlib figures)
# ---------------------------------------------------------------------------


def draw_mesh_3d(
    fig: Figure,
    vertices: np.ndarray,
    tets: np.ndarray,
    *,
    title: str = "",
    face_colors: Optional[np.ndarray] = None,
) -> None:
    fig.clear()
    ax = fig.add_subplot(111, projection="3d")
    faces = boundary_faces(tets)
    fv = face_vectors(vertices, faces)
    if face_colors is None:
        colors = [0.5, 0.5, 0.5]
        collection = mplot3d.art3d.Poly3DCollection(
            fv,
            facecolor=colors,
            lw=0.5,
            edgecolor=[0, 0, 0],
            alpha=0.85,
        )
    else:
        collection = mplot3d.art3d.Poly3DCollection(
            fv,
            facecolors=face_colors,
            lw=0.3,
            edgecolor=[0.1, 0.1, 0.1],
            alpha=0.9,
        )
    ax.add_collection3d(collection)
    scale = vertices.flatten()
    ax.auto_scale_xyz(scale, scale, scale)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()


def draw_permittivity_bars(
    fig: Figure,
    param_r: np.ndarray,
    param_i: Optional[np.ndarray] = None,
    *,
    title: str = "Permittivity by tetrahedron",
) -> None:
    fig.clear()
    if param_i is None:
        ax = fig.add_subplot(111)
        unique, counts = np.unique(param_r, return_counts=True)
        ax.bar(range(len(unique)), counts, tick_label=[f"{u:.4g}" for u in unique])
        ax.set_xlabel("Re(eps)")
        ax.set_ylabel("Tet count")
        ax.set_title(title)
    else:
        ax1 = fig.add_subplot(211)
        unique_r, counts_r = np.unique(param_r, return_counts=True)
        ax1.bar(
            range(len(unique_r)),
            counts_r,
            tick_label=[f"{u:.4g}" for u in unique_r],
            color="steelblue",
        )
        ax1.set_ylabel("Tet count")
        ax1.set_title("Re(eps) distribution")
        ax2 = fig.add_subplot(212)
        unique_i, counts_i = np.unique(param_i, return_counts=True)
        ax2.bar(
            range(len(unique_i)),
            counts_i,
            tick_label=[f"{u:.4g}" for u in unique_i],
            color="coral",
        )
        ax2.set_xlabel("Im(eps)")
        ax2.set_ylabel("Tet count")
        ax2.set_title("Im(eps) distribution")
    fig.tight_layout()


def _face_colors_by_tet_scalar(
    vertices: np.ndarray,
    tets: np.ndarray,
    tet_values: np.ndarray,
) -> np.ndarray:
    """Per-boundary-face RGBA colors from scalar per tet."""
    faces = boundary_faces(tets)
    tet_for_face = np.zeros(faces.shape[0], dtype=np.int64)
    for fi, face in enumerate(faces):
        nodes = set(face.tolist())
        for ti, tet in enumerate(tets):
            if nodes.issubset(set(tet.tolist())):
                tet_for_face[fi] = ti
                break
    vals = tet_values[tet_for_face]
    vmin, vmax = float(vals.min()), float(vals.max())
    if vmax <= vmin:
        norm = np.zeros_like(vals)
    else:
        norm = (vals - vmin) / (vmax - vmin)
    cmap = plt.cm.viridis
    return cmap(norm)


def plot_mesh(
    vertices: np.ndarray,
    faces: np.ndarray,
    *,
    save: Optional[str] = None,
    show: bool = True,
    title: Optional[str] = None,
) -> None:
    fig = plt.figure(figsize=(6, 6), frameon=False)
    ax = fig.add_subplot(111, projection="3d")
    fv = face_vectors(vertices, faces)
    ax.add_collection3d(
        mplot3d.art3d.Poly3DCollection(
            fv,
            facecolor=[0.5, 0.5, 0.5],
            lw=0.5,
            edgecolor=[0, 0, 0],
            alpha=0.8,
        )
    )
    scale = vertices.flatten()
    ax.auto_scale_xyz(scale, scale, scale)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    if title:
        ax.set_title(title)
    if save:
        plt.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Layered permittivity
# ---------------------------------------------------------------------------


def outlayer(elem: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    f = elem[:, range(4)]
    edges = np.concatenate(
        (
            f[:, [0, 1, 2]],
            f[:, [0, 1, 3]],
            f[:, [0, 2, 3]],
            f[:, [1, 2, 3]],
        ),
        axis=0,
    )
    edgesort = np.sort(edges, axis=1)
    _, uniq_idx = np.unique(edgesort, axis=0, return_index=True)
    _, inverse = np.unique(edgesort, axis=0, return_inverse=True)
    counts = np.bincount(inverse)
    qx = np.argwhere(counts == 1)
    ii = np.mod(uniq_idx[qx.ravel()], elem.shape[0])
    outind = np.unique(ii)
    return elem[outind], outind


def layers(elem: np.ndarray, depths: Sequence[int]) -> np.ndarray:
    ine = elem
    newelem = np.array([[0, 0, 0, 0, 0]], dtype=elem.dtype)
    for i in range(len(depths)):
        depth = depths[i]
        for _j in range(depth):
            oute, outi = outlayer(ine)
            ine = np.delete(ine, outi, axis=0)
            newelem = np.append(newelem, oute, axis=0)
        newelem[:, 4] = newelem[:, 4] + 1
    newelem = np.append(newelem, ine, axis=0)
    newelem = np.delete(newelem, 0, axis=0)
    return newelem


def run_layer(
    input_path: str,
    output_path: str,
    depths: Sequence[int],
    eps_r: Sequence[float],
    eps_i: Sequence[float],
) -> dict:
    if len(eps_r) != len(depths) + 1 or len(eps_i) != len(depths) + 1:
        raise ValueError("eps_r and eps_i need len(depths) + 1 values")
    data = read_h5(input_path)
    elem = np.hstack(
        [
            data["tets"].astype(np.int64) + 1,
            np.ones((data["tets"].shape[0], 1), dtype=np.int64),
        ]
    )
    elem = layers(elem, depths)
    mx = int(np.max(elem[:, 4]))
    param_r = np.zeros(elem.shape[0], dtype=np.float64)
    param_i = np.zeros(elem.shape[0], dtype=np.float64)
    for i in range(mx):
        mask = elem[:, 4] == i + 1
        param_r[mask] = eps_r[i]
        param_i[mask] = eps_i[i]
    out_path = write_h5(
        output_path,
        data["vertices"],
        elem[:, :4].astype(np.int64) - 1,
        param_r,
        param_i,
    )
    result = read_h5(out_path)
    result["n_layers"] = mx
    return result


# ---------------------------------------------------------------------------
# Gaussian random sphere
# ---------------------------------------------------------------------------


def _c_ind(n, m=None):
    if m is None:
        nn = np.floor((np.sqrt(np.asarray(n) * 8 + 1) - 1) / 2)
        nn = nn.astype(int)
        return nn, n - (nn**2 + nn) / 2
    return int(n * (n + 1) / 2 + m)


def _power_lcoef(l: np.ndarray, nu: float, lmin: int) -> np.ndarray:
    a_l = np.zeros(l.size)
    for i, value in enumerate(l):
        if value != 0:
            a_l[i] = 1.0 / value**nu
    a_l[0 : lmin + 1] = 0
    norm = np.sum(1.0 / np.arange(lmin, l.max() + 1) ** nu)
    return a_l / norm


def _mgaussian_lcoef(l: np.ndarray, ell: float, lmin: int) -> np.ndarray:
    a_l = 2.0 * (l + 1) * spherical_in(l, 1.0 / ell**2) * np.exp(-1.0 / ell**2)
    a_l[0 : lmin + 1] = 0.0
    return a_l / np.sum(a_l)


def _coef_std(a_l: np.ndarray, beta: float, l: np.ndarray, lmax: int) -> np.ndarray:
    std = beta * np.sqrt(a_l)
    cind = np.arange(lmax * (lmax + 1) // 2, dtype=int)
    ll, mm = _c_ind(cind)
    for i in cind:
        j = int(ll[i] * (ll[i] + 1) / 2)
        std[i] = std[j] * np.sqrt(
            2.0 * factorial(ll[i] - mm[i]) / factorial(ll[i] + mm[i])
        )
    return std


def _sample_gaussian_sphere_coef(
    std: np.ndarray, lmax: int, lmin: int
) -> Tuple[np.ndarray, np.ndarray]:
    a_lm = normal(size=std.size)
    b_lm = np.zeros(std.size)
    for mm in range(1, lmax):
        for ll in range(max(mm, lmin), lmax):
            i = int(_c_ind(ll, mm))
            a_lm[i] = normal() * std[i]
            b_lm[i] = normal() * std[i]
    return a_lm, b_lm


def _legendre_table(lmax: int, z: float) -> np.ndarray:
    legp = np.zeros((lmax, lmax))
    for ll in range(lmax):
        for mm in range(ll + 1):
            legp[mm, ll] = lpmv(mm, ll, z)
    return legp


def _r_gsphere(
    a_lm: np.ndarray,
    b_lm: np.ndarray,
    z: float,
    phi: float,
    beta: float,
    lmin: int,
    lmax: int,
) -> float:
    legp = _legendre_table(lmax, z)
    cphi = np.cos(np.arange(0, lmax + 1) * phi)
    sphi = np.sin(np.arange(0, lmax + 1) * phi)
    s = 0.0
    for ll in range(lmin, lmax):
        ii = int(_c_ind(ll, 0))
        s += legp[ll, 0] * a_lm[ii]
    for mm in range(1, lmax):
        for ll in range(max(mm, lmin), lmax):
            ii = int(_c_ind(ll, mm))
            s += legp[mm, ll] * (a_lm[ii] * cphi[mm] + b_lm[ii] * sphi[mm])
    return float(np.exp(s - 0.5 * beta**2))


def _deform_gsphere(
    vertices: np.ndarray,
    a_lm: np.ndarray,
    b_lm: np.ndarray,
    beta: float,
    lmin: int,
    lmax: int,
) -> np.ndarray:
    out = np.zeros_like(vertices)
    for i in range(vertices.shape[0]):
        x, y, z = vertices[i]
        phi = np.arctan2(y, x)
        norm = np.linalg.norm(vertices[i])
        z_n = z / norm if norm > 0 else 0.0
        r = _r_gsphere(a_lm, b_lm, z_n, phi, beta, lmin, lmax)
        nu = np.sqrt(max(0.0, 1.0 - z_n**2))
        out[i] = [r * nu * np.cos(phi), r * nu * np.sin(phi), r * z_n]
    return out


def run_gsphere(
    output: str,
    seed_val: int,
    refinement: float,
    append_seed: bool,
    subdivisions: int,
    sigma: float,
    nu: float,
    gamma_deg: float,
    correlation: int,
    lmin: int,
    lmax: int,
    refr: str,
) -> dict:
    gamma = np.deg2rad(gamma_deg)
    ell = 2.0 * np.sin(0.5 * gamma)
    beta = np.sqrt(np.log(sigma**2 + 1.0))
    ref_ind = complex(refr)
    eps_r = float(np.real(ref_ind**2))
    eps_i = float(np.imag(ref_ind**2))

    cind = np.arange(lmax * (lmax + 1) // 2, dtype=int)
    l_arr, _ = _c_ind(cind)
    if correlation == 1:
        a_l = _power_lcoef(l_arr, nu, lmin)
    else:
        a_l = _mgaussian_lcoef(l_arr, ell, lmin)

    seed(seed_val)
    std = _coef_std(a_l, beta, l_arr, lmax)

    base = trimesh.creation.icosphere(subdivisions=subdivisions, radius=1.0)
    a_lm, b_lm = _sample_gaussian_sphere_coef(std, lmax, lmin)
    new_vertices = _deform_gsphere(base.vertices, a_lm, b_lm, beta, lmin, lmax)
    surface = trimesh.Trimesh(vertices=new_vertices, faces=base.faces, process=True)

    nodes, tets = tetrahedralize(surface.vertices, surface.faces, refinement)
    param_r = np.full(tets.shape[0], eps_r)
    param_i = np.full(tets.shape[0], eps_i)

    out_base = output
    if append_seed:
        out_base = f"{output}{seed_val}"

    path = write_h5(out_base, nodes, tets, param_r, param_i)
    result = read_h5(path)
    result["ka_hint"] = 2 * np.pi * effective_radius(nodes, tets)
    return result


# ---------------------------------------------------------------------------
# Gaussian random ellipsoid (Muinonen & Pieniluoma 2011)
# ---------------------------------------------------------------------------


def _corr(dist: np.ndarray, ell: float) -> np.ndarray:
    r = np.asarray(dist, dtype=np.float64) / ell
    out = np.exp(-(r**2))
    out[r >= 10.0] = 0.0
    return out


def _gaussian_field(cov: np.ndarray) -> np.ndarray:
    try:
        chol = np.linalg.cholesky(cov)
        return chol @ normal(size=cov.shape[0])
    except np.linalg.LinAlgError:
        eigval, eigvec = np.linalg.eigh(cov)
        eigval = np.maximum(eigval, 0.0)
        return eigvec @ (np.sqrt(eigval) * normal(size=cov.shape[0]))


def _sphere_points(n: int) -> np.ndarray:
    idx = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * idx / n)
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    theta = 2.0 * np.pi * idx / golden
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack([x, y, z])


def _gen_ellipsoid_surface(
    a: float, b: float, c: float, n_vertices: int
) -> trimesh.Trimesh:
    pts = _sphere_points(n_vertices)
    hull = ConvexHull(pts)
    mesh = trimesh.Trimesh(vertices=pts.copy(), faces=hull.simplices, process=False)
    mesh.vertices[:, 0] *= a
    mesh.vertices[:, 1] *= b
    mesh.vertices[:, 2] *= c
    if mesh.volume < 0:
        mesh.invert()
    mesh.process(validate=True)
    return mesh


def _deform_surf(vertices: np.ndarray, ell: float, beta: float) -> np.ndarray:
    n = vertices.shape[0]
    if n > 2000:
        import warnings

        warnings.warn(
            f"deform_surf builds an {n}x{n} covariance matrix; n>2000 may be slow",
            stacklevel=2,
        )
    dist = squareform(pdist(vertices))
    cov = (beta**2) * _corr(dist, ell)
    return _gaussian_field(cov)


def _deform_ellipsoid(
    mesh: trimesh.Trimesh, heights: np.ndarray, h: float
) -> trimesh.Trimesh:
    vertices = mesh.vertices.copy()
    normals = mesh.vertex_normals
    for i in range(vertices.shape[0]):
        vertices[i] += (heights[i] - h) * normals[i]
    return trimesh.Trimesh(vertices=vertices, faces=mesh.faces, process=True)


def fix_surface(
    mesh: trimesh.Trimesh,
    detail: str = "low",
    *,
    use_pymeshfix: bool = False,
) -> trimesh.Trimesh:
    if detail not in ("low", "normal", "high"):
        raise ValueError(f"detail must be low, normal, or high, got {detail!r}")

    out = mesh.copy()
    out.update_faces(out.nondegenerate_faces())
    out.merge_vertices()
    out.process(validate=True)
    if out.volume < 0:
        out.invert()

    if use_pymeshfix or os.environ.get("SCADYN_PYMESHFIX", "") == "1":
        try:
            import pymeshfix

            meshfix = pymeshfix.MeshFix(
                np.ascontiguousarray(out.vertices, dtype=np.float64),
                np.ascontiguousarray(out.faces, dtype=np.int64),
            )
            meshfix.repair()
            out = trimesh.Trimesh(
                vertices=np.asarray(meshfix.points, dtype=np.float64),
                faces=np.asarray(meshfix.faces, dtype=np.int64),
                process=True,
            )
        except ImportError as exc:
            raise RuntimeError(
                "pymeshfix requested but not installed; "
                "pip install pymeshfix or omit --pymeshfix"
            ) from exc
        except Exception as exc:
            raise RuntimeError(f"pymeshfix repair failed: {exc}") from exc

    out.update_faces(out.nondegenerate_faces())
    out.merge_vertices()
    out.process(validate=True)
    if out.volume < 0:
        out.invert()
    return out


def _tetrahedralize_with_retry(
    surface: trimesh.Trimesh,
    refinement: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Tetrahedralize with Laplacian presmoothing (avoids TetGen native crashes)."""
    last_err: Optional[Exception] = None
    for smooth_iters in (10, 15, 20, 25):
        candidate = surface.copy()
        trimesh.smoothing.filter_laplacian(candidate, iterations=smooth_iters)
        vertices, faces = validate_surface_mesh(candidate.vertices, candidate.faces)
        for factor in (1.0, 2.0, 4.0):
            try:
                return tetrahedralize(vertices, faces, refinement * factor)
            except Exception as exc:
                last_err = exc
    raise RuntimeError(
        "Tetrahedralization failed after presmoothing and coarser retries"
    ) from last_err


def write_msh(path: str, vertices: np.ndarray, faces: np.ndarray) -> str:
    if not path.endswith(".msh"):
        path = path + ".msh"
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int64)
    with open(path, "w", encoding="utf-8") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n")
        f.write(f"{vertices.shape[0]}\n")
        for i, v in enumerate(vertices, start=1):
            f.write(f"{i} {v[0]} {v[1]} {v[2]}\n")
        f.write("$EndNodes\n")
        f.write("$Elements\n")
        f.write(f"{faces.shape[0]}\n")
        for i, tri in enumerate(faces, start=1):
            n1, n2, n3 = (int(tri[0]) + 1, int(tri[1]) + 1, int(tri[2]) + 1)
            f.write(f"{i} 2 2 0 0 {n1} {n2} {n3}\n")
        f.write("$EndElements\n")
    return path


def _truncation_hint(aeff: float, lmbda: float) -> str:
    ka = 2.0 * np.pi * aeff / lmbda
    return f"Suggested truncation check: ka = 2*pi*aeff/lambda ~ {ka:.4g}"


def run_gellip(
    output: str,
    seed_val: int,
    refinement: float,
    append_seed: bool,
    *,
    sigma: float = 0.125,
    ell: float = 0.35,
    a: float = 1.0,
    b: float = 0.8,
    c: float = 0.5,
    n_vertices: int = 1000,
    aeff: Optional[float] = 5.2,
    lmbda: Optional[float] = 0.5,
    refr: str = "1.686+0.0312j",
    detail: str = "low",
    write_msh_surface: bool = False,
    use_pymeshfix: bool = False,
) -> dict:
    if n_vertices < 4:
        raise ValueError("n_vertices must be at least 4")
    if a <= 0 or b <= 0 or c <= 0:
        raise ValueError("semiaxes a, b, c must be positive")

    h = a * (c / a) ** 2
    sigma_eff = sigma / h
    beta = np.sqrt(np.log(sigma_eff**2 + 1.0))
    ref_ind = complex(refr)
    eps_r = float(np.real(ref_ind**2))
    eps_i = float(np.imag(ref_ind**2))

    seed(seed_val)
    base_surface = _gen_ellipsoid_surface(a, b, c, n_vertices)
    roughness_field = _deform_surf(base_surface.vertices, ell, beta)

    nodes: Optional[np.ndarray] = None
    tets: Optional[np.ndarray] = None
    surface: Optional[trimesh.Trimesh] = None
    applied_amp = 1.0
    for amp in (1.0, 0.85, 0.7, 0.5, 0.35):
        heights = h + amp * roughness_field
        trial = _deform_ellipsoid(base_surface, heights, h)
        trial = fix_surface(trial, detail, use_pymeshfix=use_pymeshfix)
        try:
            nodes, tets = _tetrahedralize_with_retry(trial, refinement)
            surface = trial
            applied_amp = amp
            break
        except RuntimeError:
            continue

    if nodes is None or tets is None or surface is None:
        raise RuntimeError(
            "gellip tetrahedralization failed; try lower --sigma, -n 500, or -r 0.4"
        )

    if applied_amp < 1.0:
        import warnings

        warnings.warn(
            f"gellip: reduced roughness to {applied_amp:.2f}x for a tetrahedralizable surface",
            stacklevel=2,
        )

    scale = 1.0
    if aeff is not None:
        scale = aeff / effective_radius(nodes, tets)
        nodes = nodes * scale

    param_r = np.full(tets.shape[0], eps_r)
    param_i = np.full(tets.shape[0], eps_i)

    out_base = output
    if append_seed:
        out_base = f"{output}{seed_val}"

    path = write_h5(out_base, nodes, tets, param_r, param_i)
    result = read_h5(path)
    result["ka_hint"] = 2.0 * np.pi * effective_radius(nodes, tets)
    result["roughness_amp"] = applied_amp
    if lmbda is not None and aeff is not None:
        result["truncation_hint"] = _truncation_hint(aeff, lmbda)

    if write_msh_surface:
        msh_path = write_msh(f"{out_base}_s", surface.vertices * scale, surface.faces)
        result["surface_msh"] = msh_path

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_csv_floats(text: str) -> List[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def _parse_csv_ints(text: str) -> List[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def cmd_layer(args: argparse.Namespace) -> None:
    result = run_layer(
        args.input,
        args.output,
        _parse_csv_ints(args.depths),
        _parse_csv_floats(args.eps_r),
        _parse_csv_floats(args.eps_i),
    )
    print(f"Wrote layered mesh to {result['path']} ({result['n_layers']} layers)")


def cmd_info(args: argparse.Namespace) -> None:
    print(format_mesh_info(read_h5(args.input)))


def cmd_plot(args: argparse.Namespace) -> None:
    data = read_h5(args.input)
    faces = boundary_faces(data["tets"])
    plot_mesh(
        data["vertices"],
        faces,
        save=args.save,
        show=not args.no_show,
        title=os.path.basename(data["path"]),
    )


def cmd_gsphere(args: argparse.Namespace) -> None:
    result = run_gsphere(
        args.output,
        args.seed,
        args.refinement,
        args.append_seed,
        args.subdivisions,
        args.sigma,
        args.nu,
        args.gamma_deg,
        args.correlation,
        args.lmin,
        args.lmax,
        args.refr,
    )
    print(f"Wrote {result['path']} — {result['vertices'].shape[0]} nodes, "
          f"{result['tets'].shape[0]} tets")
    print(f"Suggested l_max check: ka ~ {result['ka_hint']}")
    if args.preview is not None:
        preview_path = args.preview if args.preview else result["path"].replace(".h5", ".png")
        plot_mesh(
            result["vertices"],
            boundary_faces(result["tets"]),
            save=preview_path,
            show=False,
        )
        print(f"Wrote preview {preview_path}")


def cmd_gellip(args: argparse.Namespace) -> None:
    result = run_gellip(
        args.output,
        args.seed,
        args.refinement,
        args.append_seed,
        sigma=args.sigma,
        ell=args.ell,
        a=args.a,
        b=args.b,
        c=args.c,
        n_vertices=args.n,
        aeff=None if args.no_aeff else args.aeff,
        lmbda=args.lmbda,
        refr=args.refr,
        detail=args.detail,
        write_msh_surface=args.write_msh,
        use_pymeshfix=args.pymeshfix,
    )
    print(f"Wrote {result['path']} — {result['vertices'].shape[0]} nodes, "
          f"{result['tets'].shape[0]} tets")
    print(f"Suggested l_max check: ka ~ {result['ka_hint']:.6g}")
    if result.get("roughness_amp", 1.0) < 1.0:
        print(
            f"Note: roughness reduced to {result['roughness_amp']:.2f}x "
            "so TetGen could mesh the surface"
        )
    if "truncation_hint" in result:
        print(result["truncation_hint"])
    if "surface_msh" in result:
        print(f"Wrote surface {result['surface_msh']}")
    if args.preview is not None:
        preview_path = args.preview if args.preview else result["path"].replace(".h5", ".png")
        plot_mesh(
            result["vertices"],
            boundary_faces(result["tets"]),
            save=preview_path,
            show=False,
        )
        print(f"Wrote preview {preview_path}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scadyn_mesh",
        description="Mesh generation and HDF5 utilities for scadyn.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    p_info = sub.add_parser("info", help="Summarize a mesh .h5 file")
    p_info.add_argument("-i", "--input", required=True)
    p_info.set_defaults(func=cmd_info)

    p_plot = sub.add_parser("plot", help="Preview a volumetric mesh")
    p_plot.add_argument("-i", "--input", required=True)
    p_plot.add_argument("--save")
    p_plot.add_argument("--no-show", action="store_true")
    p_plot.add_argument("--surface-only", action="store_true", default=True)
    p_plot.set_defaults(func=cmd_plot)

    p_layer = sub.add_parser("layer", help="Assign layered permittivity")
    p_layer.add_argument("-i", "--input", required=True)
    p_layer.add_argument("-o", "--output", required=True)
    p_layer.add_argument("--depths", default="3,8")
    p_layer.add_argument("--eps-r", default="2.128,2.128,2.0")
    p_layer.add_argument("--eps-i", default="0.003,0.003,0.003")
    p_layer.set_defaults(func=cmd_layer)

    p_gs = sub.add_parser("gsphere", help="Gaussian random sphere")
    p_gs.add_argument("-o", "--output", default="mesh")
    p_gs.add_argument("-i", "--seed", type=int, default=1)
    p_gs.add_argument("-r", "--refinement", type=float, default=0.2)
    p_gs.add_argument("--append-seed", action="store_true")
    p_gs.add_argument("--subdivisions", type=int, default=3)
    p_gs.add_argument("--sigma", type=float, default=0.285)
    p_gs.add_argument("--nu", type=float, default=4.5)
    p_gs.add_argument("--gamma-deg", type=float, default=35.0)
    p_gs.add_argument("--correlation", type=int, choices=(1, 2), default=1)
    p_gs.add_argument("--lmin", type=int, default=2)
    p_gs.add_argument("--lmax", type=int, default=10)
    p_gs.add_argument("--refr", default="1.686+0.0312j")
    p_gs.add_argument("--preview", nargs="?", const=True)
    p_gs.set_defaults(func=cmd_gsphere)

    p_ge = sub.add_parser("gellip", help="Gaussian random ellipsoid")
    p_ge.add_argument("-o", "--output", default="mesh")
    p_ge.add_argument("-i", "--seed", type=int, default=1)
    p_ge.add_argument("-r", "--refinement", type=float, default=0.2)
    p_ge.add_argument("--append-seed", action="store_true")
    p_ge.add_argument("--sigma", type=float, default=0.125)
    p_ge.add_argument("--ell", type=float, default=0.35, help="Correlation length")
    p_ge.add_argument("--a", type=float, default=1.0, help="Semiaxis a")
    p_ge.add_argument("--b", type=float, default=0.8, help="Semiaxis b")
    p_ge.add_argument("--c", type=float, default=0.5, help="Semiaxis c")
    p_ge.add_argument("-n", type=int, default=1000, help="Surface vertex count")
    p_ge.add_argument("--aeff", type=float, default=5.2, help="Target effective radius")
    p_ge.add_argument("--no-aeff", action="store_true", help="Skip volume scaling to a_eff")
    p_ge.add_argument("--lmbda", type=float, default=0.5, help="Wavelength for ka hint")
    p_ge.add_argument("--refr", default="1.686+0.0312j")
    p_ge.add_argument(
        "--detail",
        choices=("low", "normal", "high"),
        default="low",
        help="Surface repair level",
    )
    p_ge.add_argument("--write-msh", action="store_true", help="Write scaled surface .msh")
    p_ge.add_argument(
        "--pymeshfix",
        action="store_true",
        help="Use pymeshfix repair (can crash on some platforms; off by default)",
    )
    p_ge.add_argument("--preview", nargs="?", const=True)
    p_ge.set_defaults(func=cmd_gellip)

    return p


def _mesh_script_path() -> str:
    return os.path.abspath(__file__)


def _mesh_output_base(output: str, seed_val: int, append_seed: bool) -> str:
    return f"{output}{seed_val}" if append_seed else output


def _enrich_mesh_result(data: dict) -> dict:
    data["ka_hint"] = 2.0 * np.pi * effective_radius(data["vertices"], data["tets"])
    return data


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------


class PlotPanel(ttk.Frame):
    """Embedded matplotlib figure with toolbar."""

    def __init__(self, master: tk.Misc, **kwargs) -> None:
        super().__init__(master, **kwargs)
        self.figure = Figure(figsize=(6, 5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def refresh(self) -> None:
        self.canvas.draw_idle()


class ScadynMeshGUI(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("scadyn mesh tools")
        self.minsize(960, 640)
        self._busy = False
        self._current_mesh: Optional[dict] = None

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        self._build_view_tab()
        self._build_info_tab()
        self._build_gsphere_tab()
        self._build_layer_tab()
        self._build_gellip_tab()

        self.status = ttk.Label(self, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(fill=tk.X, side=tk.BOTTOM, padx=6, pady=(0, 6))

    def _set_status(self, text: str) -> None:
        self.status.config(text=text)

    def _run_bg(self, work: Callable[[], Any], on_ok: Callable[[Any], None]) -> None:
        if self._busy:
            messagebox.showwarning("Busy", "Wait for the current task to finish.")
            return
        self._busy = True
        self._set_status("Working…")

        def worker() -> None:
            err: Optional[Exception] = None
            result = None
            try:
                result = work()
            except Exception as exc:
                err = exc
            self.after(0, lambda: self._bg_done(err, result, on_ok))

        threading.Thread(target=worker, daemon=True).start()

    def _run_cli_subprocess(
        self,
        cli_args: Sequence[str],
        output_base: str,
        on_ok: Callable[[dict], None],
    ) -> None:
        """Run mesh CLI in a child process (VTK/TetGen are not thread-safe)."""
        if self._busy:
            messagebox.showwarning("Busy", "Wait for the current task to finish.")
            return
        self._busy = True
        self._set_status("Working…")
        script = _mesh_script_path()

        def worker() -> Tuple[Optional[dict], Optional[Exception]]:
            cmd = [sys.executable, script, *cli_args]
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
                if proc.returncode != 0:
                    msg = (proc.stderr or proc.stdout or "").strip()
                    if not msg:
                        msg = f"mesh command failed with exit code {proc.returncode}"
                    return None, RuntimeError(msg)
                data = _enrich_mesh_result(read_h5(output_base))
                data["cli_log"] = proc.stdout
                return data, None
            except Exception as exc:
                return None, exc

        def thread_fn() -> None:
            result, err = worker()
            if err is not None:
                self.after(0, lambda e=err: self._bg_done(e, None, on_ok))
            else:
                self.after(0, lambda r=result: self._bg_done(None, r, on_ok))

        threading.Thread(target=thread_fn, daemon=True).start()

    def _bg_done(
        self,
        err: Optional[Exception],
        result: Any,
        on_ok: Callable[[Any], None],
    ) -> None:
        self._busy = False
        if err is not None:
            self._set_status("Error")
            messagebox.showerror("Error", str(err))
            return
        on_ok(result)
        self._set_status("Done")

    @staticmethod
    def _labeled_entry(
        parent: ttk.Frame, label: str, default: str, row: int
    ) -> tk.StringVar:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W, padx=4, pady=3)
        var = tk.StringVar(value=default)
        ttk.Entry(parent, textvariable=var, width=36).grid(
            row=row, column=1, sticky=tk.EW, padx=4, pady=3
        )
        parent.columnconfigure(1, weight=1)
        return var

    def _browse_file(
        self, var: tk.StringVar, *, save: bool = False, ext: str = ".h5"
    ) -> None:
        if save:
            path = filedialog.asksaveasfilename(
                defaultextension=ext,
                filetypes=[("HDF5 mesh", "*.h5"), ("All", "*.*")],
            )
        else:
            path = filedialog.askopenfilename(
                filetypes=[("HDF5 mesh", "*.h5"), ("All", "*.*")]
            )
        if path:
            if path.endswith(".h5"):
                var.set(path[:-3])
            else:
                var.set(path)

    def _show_mesh_3d(
        self,
        panel: PlotPanel,
        data: dict,
        *,
        title: str = "",
        color_by_eps: bool = False,
    ) -> None:
        colors = None
        if color_by_eps and "param_r" in data:
            colors = _face_colors_by_tet_scalar(
                data["vertices"], data["tets"], data["param_r"]
            )
        draw_mesh_3d(panel.figure, data["vertices"], data["tets"], title=title, face_colors=colors)
        panel.refresh()

    def _show_eps_chart(self, panel: PlotPanel, data: dict) -> None:
        if "param_r" not in data:
            panel.figure.clear()
            ax = panel.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No param_r in mesh", ha="center", va="center")
            panel.refresh()
            return
        draw_permittivity_bars(
            panel.figure,
            data["param_r"],
            data.get("param_i"),
            title=os.path.basename(data["path"]),
        )
        panel.refresh()

    # --- View tab ---

    def _build_view_tab(self) -> None:
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="View mesh")

        paned = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(paned, width=320)
        right = ttk.Frame(paned)
        paned.add(left, weight=0)
        paned.add(right, weight=1)

        opts = ttk.LabelFrame(left, text="Options", padding=8)
        opts.pack(fill=tk.BOTH, expand=True)

        self.view_mesh = self._labeled_entry(opts, "Mesh file (.h5)", "mesh", 0)
        brow = ttk.Frame(opts)
        brow.grid(row=1, column=0, columnspan=2, sticky=tk.EW, pady=4)
        ttk.Button(brow, text="Browse…", command=lambda: self._browse_file(self.view_mesh)).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(brow, text="Plot 3D", command=self._on_view_plot).pack(side=tk.LEFT, padx=2)
        ttk.Button(brow, text="Permittivity chart", command=self._on_view_eps).pack(
            side=tk.LEFT, padx=2
        )

        self.view_color_eps = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            opts,
            text="Color surface by Re(eps)",
            variable=self.view_color_eps,
        ).grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=4)

        self.view_plot = PlotPanel(right)
        self.view_plot.pack(fill=tk.BOTH, expand=True)

    def _load_view_mesh(self) -> dict:
        data = read_h5(self.view_mesh.get().strip())
        self._current_mesh = data
        return data

    def _on_view_plot(self) -> None:
        def work():
            return self._load_view_mesh()

        def ok(data):
            self._show_mesh_3d(
                self.view_plot,
                data,
                title=os.path.basename(data["path"]),
                color_by_eps=self.view_color_eps.get(),
            )

        self._run_bg(work, ok)

    def _on_view_eps(self) -> None:
        def work():
            return self._load_view_mesh()

        def ok(data):
            self._show_eps_chart(self.view_plot, data)

        self._run_bg(work, ok)

    # --- Info tab ---

    def _build_info_tab(self) -> None:
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Mesh info")

        top = ttk.Frame(tab, padding=8)
        top.pack(fill=tk.X)
        ttk.Label(top, text="Mesh file").grid(row=0, column=0, sticky=tk.W, padx=4)
        self.info_mesh = tk.StringVar(value="mesh")
        ttk.Entry(top, textvariable=self.info_mesh, width=40).grid(
            row=0, column=1, sticky=tk.EW, padx=4
        )
        top.columnconfigure(1, weight=1)
        ttk.Button(top, text="Browse…", command=lambda: self._browse_file(self.info_mesh)).grid(
            row=0, column=2, padx=4
        )
        ttk.Button(top, text="Show info", command=self._on_info).grid(row=0, column=3, padx=4)

        self.info_text = scrolledtext.ScrolledText(tab, height=20, font=("Consolas", 10))
        self.info_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

    def _on_info(self) -> None:
        def work():
            return format_mesh_info(read_h5(self.info_mesh.get().strip()))

        def ok(text):
            self.info_text.delete("1.0", tk.END)
            self.info_text.insert(tk.END, text)

        self._run_bg(work, ok)

    # --- G-sphere tab ---

    def _build_gsphere_tab(self) -> None:
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="G-sphere")

        paned = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        left_outer = ttk.Frame(paned, width=340)
        right = ttk.Frame(paned)
        paned.add(left_outer, weight=0)
        paned.add(right, weight=1)

        canvas = tk.Canvas(left_outer, highlightthickness=0)
        scroll = ttk.Scrollbar(left_outer, orient=tk.VERTICAL, command=canvas.yview)
        left = ttk.Frame(canvas)
        left.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=left, anchor=tk.NW)
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        opts = ttk.LabelFrame(left, text="Parameters", padding=8)
        opts.pack(fill=tk.X, padx=4, pady=4)

        r = 0
        self.gs_output = self._labeled_entry(opts, "Output base name", "mesh", r)
        r += 1
        self.gs_seed = self._labeled_entry(opts, "RNG seed", "1", r)
        r += 1
        self.gs_refine = self._labeled_entry(opts, "Refinement (max tet vol)", "0.2", r)
        r += 1
        self.gs_append = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts, text="Append seed to output name", variable=self.gs_append).grid(
            row=r, column=0, columnspan=2, sticky=tk.W, padx=4, pady=2
        )
        r += 1
        self.gs_subdiv = self._labeled_entry(opts, "Icosphere subdivisions", "3", r)
        r += 1
        self.gs_sigma = self._labeled_entry(opts, "Sigma (radial std)", "0.285", r)
        r += 1
        self.gs_nu = self._labeled_entry(opts, "Nu (power-law index)", "4.5", r)
        r += 1
        self.gs_gamma = self._labeled_entry(opts, "Gamma (deg, Gauss corr)", "35.0", r)
        r += 1
        self.gs_corr = self._labeled_entry(opts, "Correlation (1=power, 2=Gauss)", "1", r)
        r += 1
        self.gs_lmin = self._labeled_entry(opts, "l_min", "2", r)
        r += 1
        self.gs_lmax = self._labeled_entry(opts, "l_max", "10", r)
        r += 1
        self.gs_refr = self._labeled_entry(opts, "Refractive index", "1.686+0.0312j", r)
        r += 1

        ttk.Button(opts, text="Generate & plot", command=self._on_gsphere).grid(
            row=r, column=0, columnspan=2, sticky=tk.EW, pady=8
        )

        plot_frame = ttk.Frame(right)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        self.gs_plot = PlotPanel(plot_frame)
        self.gs_plot.pack(fill=tk.BOTH, expand=True)
        self.gs_info = scrolledtext.ScrolledText(plot_frame, height=6, font=("Consolas", 9))
        self.gs_info.pack(fill=tk.X, padx=4, pady=4)

    def _on_gsphere(self) -> None:
        output = self.gs_output.get().strip()
        seed_val = int(self.gs_seed.get().strip())
        append = self.gs_append.get()
        out_base = _mesh_output_base(output, seed_val, append)
        cli = [
            "gsphere",
            "-o",
            output,
            "-i",
            str(seed_val),
            "-r",
            self.gs_refine.get().strip(),
            "--subdivisions",
            self.gs_subdiv.get().strip(),
            "--sigma",
            self.gs_sigma.get().strip(),
            "--nu",
            self.gs_nu.get().strip(),
            "--gamma-deg",
            self.gs_gamma.get().strip(),
            "--correlation",
            self.gs_corr.get().strip(),
            "--lmin",
            self.gs_lmin.get().strip(),
            "--lmax",
            self.gs_lmax.get().strip(),
            "--refr",
            self.gs_refr.get().strip(),
        ]
        if append:
            cli.append("--append-seed")

        def ok(data):
            self._current_mesh = data
            self.view_mesh.set(data["path"].replace(".h5", ""))
            self._show_mesh_3d(
                self.gs_plot,
                data,
                title=f"G-sphere — {os.path.basename(data['path'])}",
            )
            info = format_mesh_info(data) + f"\nka hint: {data['ka_hint']:.6g}"
            self.gs_info.delete("1.0", tk.END)
            self.gs_info.insert(tk.END, info)

        self._run_cli_subprocess(cli, out_base, ok)

    # --- Layer tab ---

    def _build_layer_tab(self) -> None:
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Layer permittivity")

        paned = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(paned, width=340)
        right = ttk.Frame(paned)
        paned.add(left, weight=0)
        paned.add(right, weight=1)

        opts = ttk.LabelFrame(left, text="Parameters", padding=8)
        opts.pack(fill=tk.BOTH, expand=True)

        self.ly_in = self._labeled_entry(opts, "Input mesh", "mesh", 0)
        self.ly_out = self._labeled_entry(opts, "Output mesh", "layermesh", 1)
        self.ly_depths = self._labeled_entry(opts, "Layer depths (comma)", "3,8", 2)
        self.ly_eps_r = self._labeled_entry(opts, "eps_r per layer", "2.128,2.128,2.0", 3)
        self.ly_eps_i = self._labeled_entry(opts, "eps_i per layer", "0.003,0.003,0.003", 4)

        brow = ttk.Frame(opts)
        brow.grid(row=5, column=0, columnspan=2, sticky=tk.EW, pady=6)
        ttk.Button(brow, text="Input…", command=lambda: self._browse_file(self.ly_in)).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(brow, text="Output…", command=lambda: self._browse_file(self.ly_out, save=True)).pack(
            side=tk.LEFT, padx=2
        )

        ttk.Button(opts, text="Apply layers & plot", command=self._on_layer).grid(
            row=6, column=0, columnspan=2, sticky=tk.EW, pady=4
        )

        self.layer_notebook = ttk.Notebook(right)
        self.layer_notebook.pack(fill=tk.BOTH, expand=True)
        tab3d = ttk.Frame(self.layer_notebook)
        tab_eps = ttk.Frame(self.layer_notebook)
        self.layer_notebook.add(tab3d, text="3D (colored by eps)")
        self.layer_notebook.add(tab_eps, text="Permittivity chart")

        self.layer_plot_3d = PlotPanel(tab3d)
        self.layer_plot_3d.pack(fill=tk.BOTH, expand=True)
        self.layer_plot_eps = PlotPanel(tab_eps)
        self.layer_plot_eps.pack(fill=tk.BOTH, expand=True)

        self.layer_info = scrolledtext.ScrolledText(right, height=5, font=("Consolas", 9))
        self.layer_info.pack(fill=tk.X, padx=4, pady=4)

    def _on_layer(self) -> None:
        def work():
            return run_layer(
                self.ly_in.get().strip(),
                self.ly_out.get().strip(),
                _parse_csv_ints(self.ly_depths.get()),
                _parse_csv_floats(self.ly_eps_r.get()),
                _parse_csv_floats(self.ly_eps_i.get()),
            )

        def ok(data):
            self._current_mesh = data
            self.view_mesh.set(data["path"].replace(".h5", ""))
            self._show_mesh_3d(
                self.layer_plot_3d,
                data,
                title=f"Layered — {os.path.basename(data['path'])}",
                color_by_eps=True,
            )
            self._show_eps_chart(self.layer_plot_eps, data)
            info = format_mesh_info(data) + f"\nLayers: {data['n_layers']}"
            self.layer_info.delete("1.0", tk.END)
            self.layer_info.insert(tk.END, info)

        self._run_bg(work, ok)

    # --- G-ellip tab ---

    def _build_gellip_tab(self) -> None:
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="G-ellipsoid")

        paned = ttk.Panedwindow(tab, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        left_outer = ttk.Frame(paned, width=340)
        right = ttk.Frame(paned)
        paned.add(left_outer, weight=0)
        paned.add(right, weight=1)

        canvas = tk.Canvas(left_outer, highlightthickness=0)
        scroll = ttk.Scrollbar(left_outer, orient=tk.VERTICAL, command=canvas.yview)
        left = ttk.Frame(canvas)
        left.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=left, anchor=tk.NW)
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        opts = ttk.LabelFrame(left, text="Parameters", padding=8)
        opts.pack(fill=tk.X, padx=4, pady=4)

        r = 0
        self.ge_output = self._labeled_entry(opts, "Output base name", "mesh", r)
        r += 1
        self.ge_seed = self._labeled_entry(opts, "RNG seed", "1", r)
        r += 1
        self.ge_refine = self._labeled_entry(opts, "Refinement (max tet vol)", "0.2", r)
        r += 1
        self.ge_append = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts, text="Append seed to output name", variable=self.ge_append).grid(
            row=r, column=0, columnspan=2, sticky=tk.W, padx=4, pady=2
        )
        r += 1
        self.ge_sigma = self._labeled_entry(opts, "Sigma (roughness)", "0.125", r)
        r += 1
        self.ge_ell = self._labeled_entry(opts, "Ell (correlation length)", "0.35", r)
        r += 1
        self.ge_a = self._labeled_entry(opts, "Semiaxis a", "1.0", r)
        r += 1
        self.ge_b = self._labeled_entry(opts, "Semiaxis b", "0.8", r)
        r += 1
        self.ge_c = self._labeled_entry(opts, "Semiaxis c", "0.5", r)
        r += 1
        self.ge_n = self._labeled_entry(opts, "Surface vertices n", "1000", r)
        r += 1
        self.ge_aeff = self._labeled_entry(opts, "a_eff (volume scale)", "5.2", r)
        r += 1
        self.ge_lmbda = self._labeled_entry(opts, "Lambda (ka hint)", "0.5", r)
        r += 1
        self.ge_refr = self._labeled_entry(opts, "Refractive index", "1.686+0.0312j", r)
        r += 1
        self.ge_detail = self._labeled_entry(opts, "Detail (low/normal/high)", "low", r)
        r += 1
        self.ge_write_msh = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts, text="Write surface .msh", variable=self.ge_write_msh).grid(
            row=r, column=0, columnspan=2, sticky=tk.W, padx=4, pady=2
        )
        r += 1

        ttk.Button(opts, text="Generate & plot", command=self._on_gellip).grid(
            row=r, column=0, columnspan=2, sticky=tk.EW, pady=8
        )

        plot_frame = ttk.Frame(right)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        self.ge_plot = PlotPanel(plot_frame)
        self.ge_plot.pack(fill=tk.BOTH, expand=True)
        self.ge_info = scrolledtext.ScrolledText(plot_frame, height=6, font=("Consolas", 9))
        self.ge_info.pack(fill=tk.X, padx=4, pady=4)

    def _on_gellip(self) -> None:
        output = self.ge_output.get().strip()
        seed_val = int(self.ge_seed.get().strip())
        append = self.ge_append.get()
        out_base = _mesh_output_base(output, seed_val, append)
        aeff_text = self.ge_aeff.get().strip()
        cli = [
            "gellip",
            "-o",
            output,
            "-i",
            str(seed_val),
            "-r",
            self.ge_refine.get().strip(),
            "--sigma",
            self.ge_sigma.get().strip(),
            "--ell",
            self.ge_ell.get().strip(),
            "--a",
            self.ge_a.get().strip(),
            "--b",
            self.ge_b.get().strip(),
            "--c",
            self.ge_c.get().strip(),
            "-n",
            self.ge_n.get().strip(),
            "--lmbda",
            self.ge_lmbda.get().strip(),
            "--refr",
            self.ge_refr.get().strip(),
            "--detail",
            self.ge_detail.get().strip(),
        ]
        if aeff_text:
            cli.extend(["--aeff", aeff_text])
        else:
            cli.append("--no-aeff")
        if append:
            cli.append("--append-seed")
        if self.ge_write_msh.get():
            cli.append("--write-msh")

        def ok(data):
            self._current_mesh = data
            self.view_mesh.set(data["path"].replace(".h5", ""))
            self._show_mesh_3d(
                self.ge_plot,
                data,
                title=f"G-ellipsoid — {os.path.basename(data['path'])}",
            )
            info = format_mesh_info(data) + f"\nka hint: {data['ka_hint']:.6g}"
            if data.get("cli_log"):
                info += "\n" + data["cli_log"].strip()
            if aeff_text and self.ge_lmbda.get().strip():
                info += "\n" + _truncation_hint(float(aeff_text), float(self.ge_lmbda.get()))
            msh_path = f"{out_base}_s.msh"
            if self.ge_write_msh.get() and os.path.isfile(msh_path):
                info += f"\nSurface: {msh_path}"
            self.ge_info.delete("1.0", tk.END)
            self.ge_info.insert(tk.END, info)

        self._run_cli_subprocess(cli, out_base, ok)


def launch_gui() -> None:
    app = ScadynMeshGUI()
    app.mainloop()


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] in ("--gui", "-g"):
        launch_gui()
    else:
        main()
