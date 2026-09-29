
"""
Single-Shell Free Water Tensor Model Fitting

This script fits the single-shell free water DTI model (SS-FW-DTI) on
preprocessed diffusion-weighted imaging data using a custom FreewaterRunner
module, and saves the resulting parametric maps.

Workflow:
1. Load the preprocessed DWI data, b-values, b-vectors, and brain mask
2. Crop the data to the brain bounding box and apply the mask
3. Fit the single-shell free water diffusion tensor model
4. Restore the parameter maps to the original image space
5. Save the FA, MD, AD, RD, and FW maps along with the fitting loss curve

Usage:
    python ss_fwdti_fitting.py <data_path> <dwi_file> <brainmask_file> <bval> <bvec>

Arguments:
    data_path      - directory containing the preprocessed DWI data
    dwi_file       - preprocessed DWI NIfTI file
    brainmask_file - skull-stripped brain mask file
    bval           - FSL-format bval file
    bvec           - rotated FSL-format bvec file

Dependencies: numpy, nibabel, matplotlib, dipy, pymods
"""

import os
import sys
import time

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.segment.mask import applymask, bounding_box, crop


def parse_args():
    """Parse and validate command-line arguments."""
    if len(sys.argv) < 6:
        print("Fits the single-shell FW-DTI model on single-shell diffusion data.")
        print("Usage:", sys.argv[0],
              "<data_path> <dwi_file> <brainmask_file> <bval> <bvec>")
        print("""data_path - data path with preprocessed dwi data
dwi_file - preprocessed dwi nifti file
brainmask_file - skull stripped brain mask file
bval - fsl format bval file
bvec - rotated fsl format bvec file""")
        sys.exit(0)

    return {
        'data_path': sys.argv[1],
        'dwi_file': sys.argv[2],
        'brainmask_file': sys.argv[3],
        'bval': sys.argv[4],
        'bvec': sys.argv[5],
    }


def load_data(data_path, dwi_file, brainmask_file, bval, bvec):
    """Load DWI data, gradient table, and brain mask from disk."""
    print('Loading DWI files')

    img = nib.load(os.path.join(data_path, dwi_file))
    data = img.get_fdata()
    affine = img.affine

    bvals, bvecs = read_bvals_bvecs(os.path.join(data_path, bval),
                                    os.path.join(data_path, bvec))
    gtab = gradient_table(bvals, bvecs=bvecs)

    mask_img = nib.load(os.path.join(data_path, brainmask_file))
    mask = mask_img.get_fdata().astype(bool)

    return data, affine, gtab, mask


def fit_model(masked_data, gtab, out_prefix):
    """Fit the single-shell free water tensor model and plot the loss curve."""
    module_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),
                               'pymods')
    if module_path not in sys.path:
        sys.path.append(module_path)
    from pymods.freewater_runner import FreewaterRunner

    print('Fitting single shell free-water diffusion tensor model')
    fw_runner = FreewaterRunner(masked_data, gtab)
    fw_runner.LOG = True  # turn on logging for this run
    fw_runner.run_model(num_iter=500, dt=0.001)
    fw_runner.plot_loss()
    plt.savefig(out_prefix + 'loss.png')

    return fw_runner


def restore_and_save_maps(fw_runner, data_shape, mins, maxs, affine, save_dir):
    """Restore parameter maps to the original image space and save as NIfTI."""
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    def to_full_volume(voxel_data):
        volume = np.full(data_shape[:3], np.nan)
        volume[mins[0]:maxs[0], mins[1]:maxs[1], mins[2]:maxs[2]] = voxel_data
        return volume

    output_maps = {
        'ssfweFA': fw_runner.get_fw_fa(),
        'ssfweMD': fw_runner.get_fw_md(),
        'ssfweAD': fw_runner.get_fw_ad(),
        'ssfweRD': fw_runner.get_fw_rd(),
        'FW': fw_runner.get_fw_map(),
    }

    for name, voxel_data in output_maps.items():
        out_img = nib.Nifti1Image(to_full_volume(voxel_data), affine)
        nib.save(out_img, os.path.join(save_dir, f'{name}.nii.gz'))


def main():
    """Main entry point for single-shell FW-DTI fitting."""
    plt.switch_backend('Agg')  # render figures without a display window
    start_time = time.time()

    args = parse_args()
    out_prefix = os.path.join(args['data_path'], 'ssfwe_')

    data, affine, gtab, mask = load_data(**args)

    # Crop to the brain bounding box and apply the mask
    mins, maxs = bounding_box(mask)
    mask_boolean = crop(mask, mins, maxs)
    cropped_volume = crop(data, mins, maxs)
    masked_data = applymask(cropped_volume, mask_boolean)

    fw_runner = fit_model(masked_data, gtab, out_prefix)

    save_dir = os.path.join(args['data_path'], 'ssFWDTI')
    restore_and_save_maps(fw_runner, data.shape, mins, maxs, affine, save_dir)

    print('Elapsed time:', time.time() - start_time, 'seconds')
    print('Done.')


if __name__ == '__main__':
    main()