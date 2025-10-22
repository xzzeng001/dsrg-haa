import h5py
import numpy as np
from textwrap import indent

def preview_data(dset, max_items=10):
    """Return a short preview string without loading huge arrays."""
    try:
        if dset.size == 0:
            return "[]"
        # Flatten and take up to max_items items
        arr = dset[()]  # load dataset (OK for modest sizes)
        if isinstance(arr, np.ndarray):
            flat = arr.ravel()
            head = flat[:max_items]
            return np.array2string(head, threshold=max_items, edgeitems=max_items)
        return repr(arr)
    except Exception as e:
        return f"<preview error: {e}>"

def print_attrs(obj, prefix="  "):
    if len(obj.attrs):
        print(prefix + "attrs:")
        for k, v in obj.attrs.items():
            print(prefix + f"  - {k} = {v}")

def walk_h5(filename, show_preview=True, max_preview_items=10):
    with h5py.File(filename, "r") as f:
        print(f"# File: {filename}")
        print_attrs(f, prefix="")

        def visitor(name, obj):
            if isinstance(obj, h5py.Group):
                print(f"\n[Group] /{name}")
                print_attrs(obj)
            elif isinstance(obj, h5py.Dataset):
                print(f"\n[Dataset] /{name}")
                print(f"  shape={obj.shape}, dtype={obj.dtype}")
                print_attrs(obj)
                if show_preview:
                    print("  preview:", preview_data(obj, max_preview_items))

        f.visititems(visitor)

if __name__ == "__main__":
    walk_h5("autocas_project.results_state.0.h5", show_preview=True, max_preview_items=10)

