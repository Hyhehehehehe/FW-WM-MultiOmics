"""
Free Water Tensor Model Fitting

This script fits the free water bi-tensor model (FW-DTI) for each subject
using DIPY, and saves the resulting parameter maps.

Workflow:
1. Load diffusion-weighted imaging data, b-values, b-vectors, and brain mask
2. Fit the free water bi-tensor model
3. Extract and post-process the FW fraction, FA, and MD maps
4. Save the output maps for each subject

Dependencies: numpy, nibabel, dipy
"""

import sys
import os
import numpy as np
import dipy.reconst.fwdti as fwdti
import nibabel as nib
from dipy.core.gradients import gradient_table

def check_output_files(subject_path, subject_id):
    """Strictly verify the existence of the four core output files"""
    required_outputs = [
        f"{subject_id}_original_fwFA.nii.gz",
        f"{subject_id}_processed_fwFA.nii.gz",
        f"{subject_id}_fw_volume.nii.gz",
        f"{subject_id}_fwMD.nii.gz"
    ]
    return all(os.path.exists(os.path.join(subject_path, f)) for f in required_outputs)

def process_subject(subject_id, error_log_file):
    """Main per-subject processing logic with a smart skipping mechanism"""
    base_path = os.path.dirname(error_log_file)
    subject_path = os.path.join(base_path, subject_id)

    # Pre-check: skip processing if complete results already exist
    if check_output_files(subject_path, subject_id):
        print(f"Subject {subject_id} results already exist, skipping")
        return

    try:
        # Dynamically build input file paths
        file_mapping = {
            'dwi': 'data.nii.gz',
            'bval': 'bvals',
            'bvec': 'data.eddy_rotated_bvecs',
            'mask': 'brain_mask.nii.gz'
        }

        # Load the diffusion-weighted imaging data
        dwi_img = nib.load(os.path.join(subject_path, file_mapping['dwi']))
        data = dwi_img.get_fdata()
        affine = dwi_img.affine

        # Build the gradient table
        bvals = np.loadtxt(os.path.join(subject_path, file_mapping['bval']))
        bvecs = np.loadtxt(os.path.join(subject_path, file_mapping['bvec']))
        gtab = gradient_table(bvals=bvals, bvecs=bvecs)

        # Load the brain mask
        mask_img = nib.load(os.path.join(subject_path, file_mapping['mask']))
        mask = mask_img.get_fdata().astype(bool)

        # Fit the free water bi-tensor model
        fwdtimodel = fwdti.FreeWaterTensorModel(gtab)
        fwdtifit = fwdtimodel.fit(data, mask=mask)

        # Parameter extraction and post-processing
        fw_volume = fwdtifit.f
        fwFA = fwdtifit.fa.copy()
        fwFA[fw_volume > 0.7] = 0  # Apply the free water fraction threshold correction

        # Save the results
        output_maps = {
            "original_fwFA": fwdtifit.fa,
            "processed_fwFA": fwFA,
            "fw_volume": fw_volume,
            "fwMD": fwdtifit.md
        }

        for name, array in output_maps.items():
            output_path = os.path.join(subject_path, f"{subject_id}_{name}.nii.gz")
            nib.save(nib.Nifti1Image(array, affine), output_path)

        # Report successful processing
        print(f"Subject {subject_id} processed successfully")

    except Exception as e:
        error_msg = f"Processing failed for {subject_id}: {str(e)}\n"
        with open(error_log_file, 'a') as f:
            f.write(error_msg)
        print(error_msg)

def main():
    """Main entry point for command-line argument handling"""
    if len(sys.argv) < 2:
        print("A directory number (e.g., 10) is required. Terminating.")
        sys.exit(1)

    a_number = sys.argv[1]
    # Base directory containing the unpacked subject folders (A_<number>).
    # Pass it as the second command-line argument; defaults to ./A_<number>.
    base_path = sys.argv[2] if len(sys.argv) > 2 else f"./A_{a_number}"
    error_log = os.path.join(base_path, "error_log.txt")

    # Dynamically obtain the subject list
    subjects = [d for d in os.listdir(base_path)
               if os.path.isdir(os.path.join(base_path, d))]

    print(f"Processing directory A_{a_number}, found {len(subjects)} subjects", flush=True)
    for subject in subjects:
        process_subject(subject, error_log)

if __name__ == "__main__":
    main()
