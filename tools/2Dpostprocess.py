import h5py
from pathlib import Path

def process_file(filename):
    fs = h5py.File(filename, 'r')

    # Extract the base filename without extension
    base_name = filename.stem

    # Enumerate the step groups actually present rather than reconstructing
    # their zero-padded name from NumberOfSteps (the padding width is fixed
    # by the exporter, not derived from the step count).
    step_names = sorted(
        name for name in fs.keys()
        if name.startswith(f'{base_name}_') and f'{name}/VTKHDF' in fs
    )

    for name in step_names:
        output_file = f'results2D/postproc/{name}.hdf'
        with h5py.File(output_file, 'w') as out:
            fs.copy(f'{name}/VTKHDF', out)

    fs.close()

# Loop over the raw per-slice files written by REIMS into results2D/
data_path = Path('results2D')
Path('results2D/postproc').mkdir(parents=True, exist_ok=True)
for hdf_file in data_path.glob('*.hdf'):
    process_file(hdf_file)
