import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os
from pathlib import Path

def baseline_als(y, lam, lam1, p, niter=10):
    """
    Automated baseline correction using asymmetric least squares smoothing.
    This is the same algorithm from the original GUI program.
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    D1 = sparse.diags([1, -1], [0, -1], shape=(L, L - 1))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W.dot(W.transpose()) + lam1 * D1.dot(D1.transpose()) + lam * D.dot(D.transpose())
        b = (W.dot(W.transpose()) + lam1 * D1.dot(D1.transpose())).dot(y)
        z = spsolve(Z, b)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

def process_spectra(input_folder, output_folder, lam=300, lam1=0.001, p=0.01):
    """
    Process all .asc files in the input folder and save baseline-corrected spectra to the output folder.
    """
    # Create output directory if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Get list of all .asc files in input directory
    asc_files = [f for f in os.listdir(input_folder) if f.endswith('.asc')]
    
    print(f"Found {len(asc_files)} .asc files to process")
    
    # Process each file
    for filename in asc_files:
        print(f"Processing {filename}...")
        
        # Load the spectrum
        input_path = os.path.join(input_folder, filename)
        data = np.loadtxt(input_path)
        
        # Extract wavelength and intensity
        wavelength = data[:, 0]
        intensity = data[:, 1]
        
        # Perform baseline correction
        baseline = baseline_als(intensity, lam, lam1, p)
        corrected_intensity = intensity - baseline
        
        # Save the corrected spectrum
        output_path = os.path.join(output_folder, filename)
        with open(output_path, 'w') as f:
            for w, i in zip(wavelength, corrected_intensity):
                f.write(f"{w}\t{i}\n")
        
        print(f"Saved corrected spectrum to {output_path}")

if __name__ == "__main__":
    # Define your input and output folders
    input_folder = r"C:\Users\ASUS\Desktop\Research\Results\Raman\05.12.2024 - Fluoride series with and without heat fold\Fluoride series without heatfold"
    output_folder = r"C:\Users\ASUS\Desktop\Research\Results\Raman\05.12.2024 - Fluoride series with and without heat fold\Fluoride series without heatfold baseline corrected"
    
    # Process all spectra with the specified parameters
    process_spectra(input_folder, output_folder, lam=300, lam1=0.001, p=0.01)
    
    print("Baseline correction completed for all spectra!")
