import ctypes
import numpy as np

he_prop_dll = ctypes.CDLL("he_prop.dll", winmode=0)   # Load the DLL

# Prototype: type, size, num_cols, table_pointer
he_prop_dll.he_prop_all.restype = None
he_prop_dll.he_prop_all.argtypes = [ ctypes.c_int, ctypes.c_int, ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="F_CONTIGUOUS")]

def he_generic(func_type, inputs, num_outputs):
    """ Generic wrapper for Fortran helium routines. `inputs` is a tuple of input arrays/scalars."""
    # Broadcast inputs so scalars and arrays map to the same shape
    b_inputs = np.broadcast_arrays(*inputs)
    shape = b_inputs[0].shape
    size  = b_inputs[0].size
    num_inputs = len(b_inputs)
    num_cols = num_inputs + num_outputs

    # Create Fortran-contiguous array to prevent transposition errors
    data = np.empty((size, num_cols), dtype=np.float64, order="F")
    
    # Fill inputs (columns 0 to num_inputs-1 in Python)
    for i in range(num_inputs): data[:, i] = np.asarray(b_inputs[i], dtype=np.float64).flatten()
        
    # Call Fortran
    he_prop_dll.he_prop_all(func_type, size, num_cols, data)

    # Extract outputs (columns num_inputs to end)
    outputs = []
    for i in range(num_outputs):
        out_col = data[:, num_inputs + i].reshape(shape)
        outputs.append(out_col)

    # Return single array if 1 output, else tuple
    if num_outputs == 1: return outputs[0]
    return tuple(outputs)

# --- 1-Output Functions ---
def r_roT(ro, t): return he_generic(1, (ro, t), num_outputs=1)
def ro_pT(p, t):  return he_generic(2, (p, t),  num_outputs=1)
def T_roP(ro, p): return he_generic(3, (ro, p), num_outputs=1)
def T_roE(ro, e): return he_generic(4, (ro, e), num_outputs=1)

# --- Complex/Multi-Output Functions ---
def droeint_droP(ro, p):
    """Returns droeint_dro, e, dedv, dedp"""
    return he_generic(5, (ro, p), num_outputs=4)

def state_roT(ro, t):
    """Returns r, p, e, c"""
    return he_generic(6, (ro, t), num_outputs=4)

def state_roP(ro, p):
    """Returns r, e, t, c"""
    return he_generic(7, (ro, p), num_outputs=4)

def state_roP_withR(ro, p, r):
    """Returns e, t, c"""
    return he_generic(8, (ro, p, r), num_outputs=3)

def state_roE(ro, e):
    """Returns p, t, c"""
    return he_generic(9, (ro, e), num_outputs=3)

def state_roE_withR(ro, e, r):
    """Returns p, t, c"""
    return he_generic(10, (ro, e, r), num_outputs=3)

def jacobian_roT(ro, t):
    """Returns cv, dPdT, dTdP, dTdRo, dRdRo, dRdT, r, p, d2TdP_dT, d2TdP_dRo"""
    return he_generic(11, (ro, t), num_outputs=10)

def dc2_roT(ro, t):
    """Returns dc2dRo, dc2dP"""
    return he_generic(12, (ro, t), num_outputs=2)

def he_prop(ro, t):
    """Returns cv, cp, mu, lambda, dPdT"""
    return he_generic(13, (ro, t), num_outputs=5)

if __name__ == "__main__":
    # --- Example Usage ---
    p_arr = np.array([3e5, 4e5])
    t_arr = np.array([300, 200])
    
    print("Testing ro_pT:")
    ro = ro_pT(p_arr, t_arr)
    print(f"ro: {ro}\n")
    
    print("Testing state_roP_withR (Mixing scalar and vectors!):")
    # Using the ro array we just found, plus scalar r=1.0
    e, t, c = state_roP_withR(ro, p_arr, r=1.0)
    print(f"e: {e}")
    print(f"t: {t}")
    print(f"c: {c}")