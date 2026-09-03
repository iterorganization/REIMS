
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import h5py
import numpy as np

# exports only copy_hdf5 function
def copy_hdf5(source_file, target_file):
    def copy_attrs(source, target):
        for key, value in source.attrs.items(): # Ignore scale attributes
            if key not in ['DIMENSION_LIST', 'REFERENCE_LIST']:
                target.attrs[key] = value

    def copy_group(source_group, target_group, dims={}):
        copy_attrs(source_group, target_group)
        for name, item in source_group.items():
            if isinstance(item, h5py.Group):
                new_group = target_group.create_group(name)
                copy_group(item, new_group, dims)
            elif isinstance(item, h5py.Dataset):
                target_dataset = target_group.create_dataset(
                    name, shape=item.shape, dtype=item.dtype)
                copy_attrs(item, target_dataset)
                dims[item.name] = {}
                for i,dim in enumerate(item.dims):
                    if len(dim) == 0: continue
                    dims[item.name][(i,dim.label)] = {k:v.name for k,v in dim.items()}
                if not item.chunks:
                    target_dataset[:] = item[:]
                    continue
                # Copy chunked dataset as contiguous dataset by iterating over chunks
                buffer = np.empty(item.shape, dtype=item.dtype)                
                seg = np.ceil(np.array(item.shape) / np.array(item.chunks)).astype(int)
                for chunk_idx in np.ndindex(*seg):
                    slices = tuple(slice(i * cs, min((i + 1) * cs, s))
                        for i, cs, s in zip(chunk_idx, item.chunks, item.shape))
                    buffer[slices] = item[slices]
                target_dataset[:] = buffer[:]
        return dims

    with h5py.File(source_file, 'r') as src, h5py.File(target_file, 'w') as tgt:
        dims = copy_group(src, tgt)
        for ds, dim in dims.items(): # Reconstruct dimension scales
            for label, scales in dim.items():
                for scale, name in scales.items():
                    tgt[ds].dims[label[0]].attach_scale(tgt[name])
                tgt[ds].dims[label[0]].label = label[1]

if __name__ == '__main__':
    # Process command line arguments if run as a script
    from argparse import ArgumentParser
    from pathlib import Path

    par = ArgumentParser(description='Copy hdf5 files with removed chunks transforming them to continuos datasets.')
    par.add_argument('in_files',          nargs='+', type=Path, help='input file or files')
    par.add_argument('-o', '--out_files', nargs='*', type=Path, help='output file or files (optional)')
    par.add_argument('-a', '--postfix', default='_copy', help='postfix to append if output files are not provided')
    args = par.parse_args()

    if args.out_files:
        if len(args.in_files) != len(args.out_files):
            raise ValueError("Number of input and output files must match when output files are specified.")
        out_files = args.out_files
    else:
        out_files = [in_file.with_stem(in_file.stem + args.postfix) for in_file in args.in_files]

    for file_in, file_out in zip(args.in_files, out_files): copy_hdf5(file_in, file_out)
