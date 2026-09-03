from importlib import reload
import h5py
import xarray as xr
import numpy as np
import h5netcdf

var_names = dict(
    P  ='Pressure, Pa',
    rho='Density, kg/m³',
    u  ='u',
    m  ='Mass flow rate, kg/s',
    T  ='Temperature, K',
    mH ='mH',
    R  ='R'
)

def to_node_str(val,max_val,str_len):
    if 1 <= val <= max_val: return str(int(val)).zfill(str_len)
    if -max_val-1 < val < 0: return str(int(max_val+val+1)).zfill(str_len)
    if 0 < val < 1: return str(int(round(val*(max_val-1)+1))).zfill(str_len)

def DataArray_jsel(arr:xr.DataArray,cmp=[],node=[],name=[],**forward_to_sel):
    '''Jacek's select: da.jsel(cmp=['bg','hl'],node=[1,.5,-1])'''
    arr_cmp = {}
    for s in arr.coords['name'].values:
        k, v = s.rsplit(' ', 1)
        if arr_cmp.get(k, (0,0))[0] < int(v): arr_cmp[k] = int(v),len(v)
    if cmp == []: cmp = list(arr_cmp.keys())
    if type(cmp) is not list: cmp = [cmp]
    #cmp = set(cmp+list(arr_cmp.keys()))
    gen_nodes = True if node == [] else False 
    names = set(name)
    for c in cmp:
        max_val,str_len = arr_cmp.get(c,(None, None))
        if max_val == None: continue
        if gen_nodes: node = [i for i in range(1,max_val+1)]
        #print(c,max_val,str_len,node)
        for val in node:
            node_str = to_node_str(val, max_val, str_len)
            if node_str != None: names.add(c + ' ' + node_str)
    names = list(names & set(arr['name'].values))
    names.sort()
    return arr.sel(name=names,**forward_to_sel)
xr.DataArray.jsel = DataArray_jsel

def Dataset_get_arr(dataset,name) -> xr.DataArray:
    arr:xr.DataArray = dataset[name]
    arr = arr.transpose(name+'_var',name+'_name','t')
    for k in arr.coords.keys():
        if k.startswith(name):
            arr = arr.rename({k:k[len(name)+1:]})
    nodes = [arr.nodes] if type(arr.nodes) == np.int32 else arr.nodes
    arr = arr.assign_coords({
        'x'   : ('name', dataset.coords[name + '_x'  ].values),
        'cmp' : ('name', dataset.coords[name + '_cmp'].values),
        'node': ('name', sum([list(range(1,i+1)) for i in nodes],[])),
    })
    arr.attrs={}
    return arr
xr.Dataset.get_arr = Dataset_get_arr

def load_h5(h5_file_name):
    # engine = 'netcdf4' or 'h5netcdf'
    with xr.open_dataset(h5_file_name,group='reims',engine='h5netcdf') as h5:
        return h5.load()
    
def DataArray_run(data,prefix):
    data2=data.copy()
    data2['name'] = xr.apply_ufunc(np.vectorize(lambda val:prefix+' '+val),data2['name'])
    data2['cmp']  = xr.apply_ufunc(np.vectorize(lambda val:prefix+' '+val),data2['cmp'])
    return data2
xr.DataArray.run = DataArray_run

def DataArray_save_animation(arr: xr.DataArray, file_name, y_label='', line_style={},
        legend=True, legend_arg={}, grid={'visible':False}, dpi=100, fps=30):
    from matplotlib import pyplot as plt
    from matplotlib import animation

    comp_names = np.unique(arr.coords['cmp'])
    comps = [arr.coords['cmp'] == cmp for cmp in comp_names]
    t_max = arr.coords['t'].max()
    if y_label == '': y_label = var_names[str(arr.coords['var'].values)]

    fig, ax = plt.subplots(dpi=dpi)
    ax.grid(**grid)
    fig.subplots_adjust(top=0.95, bottom=0.15)
    ax.set_ylim(arr.min(), arr.max())
    ax.set_xlim(arr.coords['x'].min(), arr.coords['x'].max())
    lines = [ax.plot(
        arr.isel(t=0, name=cmp).coords['x'].values,  # x
        arr.isel(t=0, name=cmp).values,              # y
        label=comp_names[i],
        **line_style.get(comp_names[i],{})
    )[0] for i, cmp in enumerate(comps)]
    if legend: ax.legend(**legend_arg)
    ax.set_xlabel('x, m')
    ax.set_ylabel(y_label)
    
    progress_ax = fig.add_axes([0.1, 0.01, 0.8, 0.03])  # [left, bottom, width, height]
    progress_ax.set_xlim(0, 1)
    progress_ax.set_ylim(0, 1)
    progress_ax.axis('off')  # Hide the axis
    progress_bar = progress_ax.barh(0.5, 0, height=1, color='lightgrey')[0]
    time_text = progress_ax.text(0.5, 0.5, '', transform=progress_ax.transAxes, 
                                 fontsize=10, verticalalignment='center', 
                                 horizontalalignment='center', color='black')
    
    def update(frame):
        for i, cmp in enumerate(comps):
            lines[i].set_ydata(arr.isel(t=frame, name=cmp).values)
        current_time = arr.coords['t'][frame].values
        time_text.set_text(f'Time = {current_time:.2f}')
        progress_bar.set_width(current_time / t_max)
        return lines + [progress_bar, time_text]
    
    ani = animation.FuncAnimation(fig, update, frames=arr.sizes['t'], blit=True)
    ani.save(file_name, writer='ffmpeg', fps=fps)
xr.DataArray.save_animation = DataArray_save_animation


# Fast read --------------------------------------------------------------------
def fast_dset_read(dset):
    data = np.empty(dset.shape, dtype=dset.dtype)
    slices = [tuple(slice(i * cs, min((i + 1) * cs, s))
        for i, cs, s in zip(chunk_idx, dset.chunks, dset.shape))
        for chunk_idx in np.ndindex(*[-(-s//c) for s,c in zip(dset.shape, dset.chunks)])]
    for s in slices: data[s] = dset[s]
    return data

def can_fast_read(self, key):
    if not self._h5ds.chunks:                   return False
    if not isinstance(self, h5netcdf.Variable): return False
    if isinstance(key, xr.core.indexing.BasicIndexer): key_tuple = key.tuple
    elif isinstance(key, tuple):                       key_tuple = key
    else:                                       return False
    if len(key_tuple) < len(self.shape):
        key_tuple = key_tuple + (slice(None),) * (len(self.shape) - len(key_tuple))
    return all(isinstance(k, slice) and k == slice(None) for k in key_tuple)

original_getitem = h5netcdf.Variable.__getitem__
def fast_getitem(self, key):
    if can_fast_read(self, key): return fast_dset_read(self._h5ds)
    return original_getitem(self, key)
h5netcdf.Variable.__getitem__ = fast_getitem
